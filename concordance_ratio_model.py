import pandas as pd 
import yaml


# Main functions

def prepare_df_for_concordance_ratio(df, validity_duration, evaluation_start_date, evaluation_end_date, patient_col, date_col):
    """
    Prepare a DataFrame for calculating the concordance score with the ratio model by processing activity dates 
    and ensuring proper boundary handling for the evaluation period.

    Steps performed
    ---------------
    1. Filter activity dates to include only those within the evaluation period 
       (between `evaluation_start_date` and `evaluation_end_date`).
    2. Calculate a "dummy date," defined as `evaluation_start_date` minus the `validity_duration`.
    3. For each patient:
       - Retain one activity date before the `evaluation_start_date` to account for prior activities:
         - If no activity date exists before `evaluation_start_date`, or the last activity date is earlier 
           than the dummy date, use the dummy date.
         - Otherwise, use the most recent activity date before `evaluation_start_date`.
         - Note: Using the dummy date does not contribute to concordant days in the evaluation period, 
           as it only ensures continuity for calculations.
    4. Add a "test date" set to one day after the `evaluation_end_date` to handle boundary conditions 
       during concordance calculations.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing patient_col, date_col, and binary column(s) for each indicator of interest.
    validity_duration : float
        The length of the validity period (in days).
    evaluation_start_date : pandas.Timestamp
        The start date of the evaluation period.
    evaluation_end_date : pandas.Timestamp
        The end date of the evaluation period.
    patient_col : str
        Name of the column containing (masked) patient identifiers.
    date_col : str
        Name of the column containing dates on which the indicator was met.

    Returns
    -------
    df_final : pandas.DataFrame
        A DataFrame filtered and augmented according to the specified rules.    
    """
    # Formulate dummy date, which is set to the evaluation_start_date minus the validity_duration
    dummy_date = evaluation_start_date - pd.Timedelta(days=validity_duration)

    # Dataframe with most recent activitiy date before the start of the evaluation period for each patient
    last_date_before_evaluation = df[df[date_col] < evaluation_start_date].groupby(patient_col)[date_col].max().reset_index()

    # If the last test day before the evaluation falls before dummy date, replace by dummy date instead    
    last_date_before_evaluation[date_col] = last_date_before_evaluation[date_col].clip(lower=dummy_date)

    # Identify patients with no dates before evaluation_start_date and assign them the dummy date
    all_patients = pd.DataFrame({patient_col: df[patient_col].unique()})
    last_date_before_evaluation = all_patients.merge(last_date_before_evaluation, on=patient_col, how='left')
    last_date_before_evaluation[date_col] = last_date_before_evaluation[date_col].fillna(dummy_date)

    # DataFrame that only contains activities that were completed within the evaluation period
    df_evaluation_period = df[(df[date_col] >= evaluation_start_date) & (df[date_col] <= evaluation_end_date)]

    # Add extra date one day after the evaluation_end_date
    df_day_after_end_date = all_patients.copy()
    df_day_after_end_date[date_col] = evaluation_end_date + pd.Timedelta(days=1)

    # Combine DataFrames
    df_combined = pd.concat([df_evaluation_period, last_date_before_evaluation, df_day_after_end_date], ignore_index=True)
    df_combined.sort_values(by=[patient_col, date_col], inplace=True)

    # Remove duplicates 
    df_combined.drop_duplicates(inplace=True)

    return df_combined
    

def calc_ratio_of_concordant_period(df, indicator, evaluation_start_date, evaluation_length, validity_duration, patient_col, date_col):
    """
    This function calculates the ratio of covered days to total days in the evaluation period. 

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing data about clinical activities, including the (masked) patient IDs, dates,
        and activity types. The DataFrame can have multiple rows (i.e., dates) per patient.
    indicator : str
        The indicator for which the concordance score will be calculated.
    evaluation_start_date : str
        The start date of the evaluation period in 'YYYY-MM-DD' format.
    evaluation_length : int
        Duration of the evaluation period in days.
    validity_duration : float
        The validity duration in days for the given indicator.
    patient_col : str
        Name of the column containing (masked) patient identifiers.
    date_col : str
        Name of the column containing dates on which the indicator was met.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with two columns: the (masked) patient ID and concordance score for the indicator of interest
    """
    # Calculate evaluation end date
    evaluation_end_date = evaluation_start_date + pd.Timedelta(days=evaluation_length - 1) 
        
    # Filter on specific indicator and only keep relevant columns
    df_filtered = df[df[indicator] == 1][[date_col, indicator, patient_col]]

    # Prepare dataframe for calculating the ratio by only keeping the dates in the evaluation period + one test date before the evaluation start date
    df_filtered = prepare_df_for_concordance_ratio(df_filtered, validity_duration, evaluation_start_date, evaluation_end_date, patient_col, date_col)

    # Shift the date_col to get the 'prev_date' for each patient
    df_filtered['prev_date'] = df_filtered.groupby(patient_col)[date_col].shift(1)

    # Drop rows where 'prev_date' is NaT; these are the first rows of each group
    df_filtered = df_filtered.dropna(subset=['prev_date'])

    # Calculate days between two activities
    df_filtered['days_between_activities'] = (df_filtered[date_col] - df_filtered['prev_date']).dt.total_seconds() / (24 * 3600)

    # Subtract the validity duration from each time between two activities 
    df_filtered['days_not_concordant'] = df_filtered['days_between_activities'] - validity_duration

    # Ensure that 'days_not_concordant' does not go below 0 (this happens when duration between two tests is shorter than the maximum coverage duration of the test)
    df_filtered['days_not_concordant'] = df_filtered['days_not_concordant'].clip(lower=0)

    # Calculate ratio of coverage 
    df_filtered['concordance_{}'.format(indicator)] = (evaluation_length - df_filtered.groupby(patient_col)['days_not_concordant'].transform('sum')) / evaluation_length

    # Return dataframe with concordance score for each patient
    df_score = df_filtered.drop_duplicates(keep='first', subset=patient_col)[[patient_col, 'concordance_{}'.format(indicator)]]

    return df_score


def calc_concordance_with_ratio(df_, evaluation_start_date, evaluation_length, indicators, patient_col, date_col, 
                                validity_periods=None, validity_periods_file="validity_periods.yml"): 
    """ 
    This function calculates the concordance score via the ratio model, as described in [$Paper Title] ($link_to__paper).
    The concordance score of a single activity is calculated as follows: 

        concordance score = [number of concordant days in evaluation period] / [number days in evaluation period]

    The total concordance is calculated by averaging the concordance scores for all individual indicators.

    Parameters
    ----------
    df_ : pandas.DataFrame
        The DataFrame containing data about clinical activities, including (masked) patient IDs, dates,
        and activity types. The DataFrame can have multiple rows (i.e., dates) per patient.
    evaluation_start_date : str or pandas.Timestamp
        The start date of the evaluation period in 'YYYY-MM-DD' format.
    evaluation_length : int
        Duration of the evaluation period in days.
    indicators : list of str or str
        List or single string indicating the indicator(s) for which to calculate the concordance score.
        Example: ['eGFR', 'HbA1c'] or 'eGFR'.
    patient_col : str
        Name of the column containing (masked) patient identifiers.
    date_col : str
        Name of the column containing dates on which the indicator was met.
    validity_periods : dict or None, optional
        Dictionary with indicator as key and validity duration in days as value. If None, data from
        `validity_periods_file` will be used.
        Default: None.
    validity_periods_file : str, optional
        Path to the YAML file containing validity periods.
        Default: 'validity_periods.yml'.

    Returns
    -------
    pandas.DataFrame
        DataFrame with concordance scores for all requested indicators (columns) for each patient (rows).
    """
    # Copy the DataFrame to avoid modifying the original
    df = df_.copy()

    # Load and validate validity periods 
    validity_periods_lower = load_validity_periods(validity_periods, validity_periods_file)

    # Validate inputs
    df, evaluation_start_date, indicators = validate_inputs(df, evaluation_start_date, evaluation_length, indicators, patient_col, date_col, validity_periods_lower)

    # Calculate concordance scores for individual indicators
    df_scores_all = pd.DataFrame(df[patient_col].unique(), columns=[patient_col]) # dataframe with all unique patients
    
    for indicator in indicators:
        # Check if there are activities for this indicator without a date
        missing_dates = df[df[indicator] == 1][date_col].isna().sum()

        if missing_dates > 0:
            print(f"WARNING: The indicator '{indicator}' has {missing_dates} activities without a date in the column '{date_col}'.\n"
                f"Activities without a date will not be included in the concordance score calculation!")

        # Get the validity duration for this indicator
        validity_duration = validity_periods_lower[indicator.lower()]

        df_scores = calc_ratio_of_concordant_period(df, indicator, evaluation_start_date, evaluation_length, validity_duration, patient_col, date_col)
        df_scores_all = pd.merge(df_scores_all, df_scores, on=patient_col, how='outer') # outer: keep patients who did not perform any tests
    
    # Calculate average concordance score in case of more than one clinical indicator
    if len(indicators) > 1:
        cols_with_scores = df_scores_all.columns[df_scores_all.columns.str.startswith('concordance')]
        df_scores_all['concordance_total'] = df_scores_all[cols_with_scores].mean(axis=1)

    return df_scores_all



# Validation functions

def load_validity_periods(validity_periods, validity_periods_file):
    """
    Loads and validates the validity periods from a dictionary or a YAML file, and returns it as a lower-case dictionary.

    Parameters
    ----------
    validity_periods : dict or None
        Dictionary with indicator : validity duration (in days). If None, the YAML file is used.
    validity_periods_file : str
        Path to the YAML file containing validity periods.

    Returns
    -------
    dict
        Validity periods with lowercase keys for case-insensitive matching.

    Raises
    ------
    FileNotFoundError
        If the validity periods file is not found.
    RuntimeError
        If an error occurs while loading the YAML file.
    TypeError
        If `validity_periods` is not a dictionary when provided.
    """
        
    if validity_periods is None:
        try:
            with open(validity_periods_file, "r") as f:
                validity_periods = yaml.safe_load(f)
            validity_periods_lower = {key.lower(): value for key, value in validity_periods.items()}
        except FileNotFoundError:
            raise FileNotFoundError(f"The file '{validity_periods_file}' could not be found.")
        except Exception as e:
            raise RuntimeError(f"An error occurred while loading '{validity_periods_file}': {e}")
    else:
        if not isinstance(validity_periods, dict):
            raise TypeError("The `validity_periods` parameter must be a dictionary.")
        validity_periods_lower = {key.lower(): value for key, value in validity_periods.items()}
    return validity_periods_lower


def validate_inputs(df, evaluation_start_date, evaluation_length, indicators, patient_col, date_col, validity_periods_lower):
    """
    Validates the inputs for calculating concordance scores.

    This function ensures that the necessary columns are present in the DataFrame, verifies the evaluation period, 
    and checks for the presence of validity periods for the specified indicators.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing data about clinical activities, including the (masked) patient IDs, dates, and activity types.
    evaluation_start_date : str
        Start date of the evaluation period in 'YYYY-MM-DD' format.
    evaluation_length : int
        Duration (in days) of the evaluation period.
    indicators : list of str or str
        List of strings or a single string indicating the indicators to validate. 
        Example: ['eGFR', 'HbA1c'] or 'eGFR'.
    patient_col : str
        Name of the column with (masked) patient identifiers.
    date_col : str
        Name of the column with dates on which the clinical activity was performed.
    validity_periods_lower : dict
        Dictionary with indicator : validity duration in days as key/value pairs.
        Keys are case-insensitive.

    Returns
    -------
    df : pandas.DataFrame
        The validated DataFrame.
    evaluation_start_date : pandas.Timestamp
        The validated and converted evaluation start date as a datetime object.
    indicators : list
        Validated list of indicators (always returned as a list, even if input is a single string).

    Raises
    ------
    TypeError
        If the input DataFrame is not a pandas DataFrame or if the `indicators` parameter is not a list or string.
    ValueError
        If required columns are missing from the DataFrame, if the `evaluation_start_date` format is invalid, 
        or if any specified indicators do not have a specified validity period.

        
    Notes
    -----
    - Matching between indicators and indicator keys in validity_periods_lower is case-insensitive.
    """

     # Validate DataFrame type
    if not isinstance(df, pd.DataFrame):
        raise TypeError("The `df` parameter must be a pandas DataFrame.")
    
    # Check if there is any data in the dataframe
    if df.empty:
        raise ValueError("The input DataFrame is empty. Please provide a valid DataFrame.")
    
    # Convert date column to datetime
    try:
        df[date_col] = pd.to_datetime(df[date_col])
    except Exception as e:
        raise ValueError(f"Could not convert the `{date_col}` column to datetime format. Ensure it contains valid dates.") from e

    # Convert evaluation_start_date to datetime
    try:
        evaluation_start_date = pd.to_datetime(evaluation_start_date, format='%Y-%m-%d')
    except Exception as e:
        raise ValueError(f"Invalid date format for `evaluation_start_date`: {evaluation_start_date}. Expected format: 'YYYY-MM-DD'.") from e

    # Check if the evaluation length is a positive number
    if not isinstance(evaluation_length, (int, float)) or evaluation_length <= 0:
        raise ValueError("The `evaluation_length` parameter must be a positive number representing the number of days.")
    
    # Check if any indicators are given 
    if not indicators:
        raise ValueError("The `indicators` parameter cannot be empty. Provide at least one indicator.")

    # Ensure indicators is a list or string
    if isinstance(indicators, str):
        indicators = [indicators]  # Convert single string to a list
    elif not isinstance(indicators, (list)):
        raise TypeError("The `indicators` parameter must be a list (one or more indicators) or string (single indicator).")
    
    # Check required columns in the DataFrame
    required_columns = [patient_col, date_col] + indicators
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"The following required columns are missing from the DataFrame: {missing_columns}")
               
    # Check if all indicators have a validity period specified
    missing_indicators = [ind for ind in indicators if ind.lower() not in validity_periods_lower]
    if missing_indicators:
        raise ValueError(f"The following indicators do not have specified validity periods: {missing_indicators}")

    return df, evaluation_start_date, indicators