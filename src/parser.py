# parser.py

import configparser
import os

def parse_input_file(filepath):
    """
    Parses the simulation configuration file.

    Args:
        filepath (str): The path to the configuration file.

    Returns:
        dict: A dictionary containing the simulation parameters.
              Values are converted to appropriate types (int, float, bool)
              where possible, otherwise kept as strings.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Configuration file not found at: {filepath}")

    params = {}
    try:
        config = configparser.ConfigParser()
        # Use allow_no_value=True if you might have key-only options (less common)
        # Add optionxform to preserve case of keys
        config.optionxform = str

        # configparser expects sections, but we can use a default section
        # Or, read line by line manually if configparser is too strict
        # Let's write a simple manual parser for more flexibility with comments and format

        with open(filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue # Skip empty lines and comments

                if '=' not in line:
                    print(f"Warning: Skipping malformed line {line_num} in {filepath}: '{line}' (no '=')")
                    continue

                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()

                # Attempt to convert value to appropriate type
                try:
                    # Try integer
                    params[key] = int(value)
                except ValueError:
                    try:
                        # Try float
                        params[key] = float(value)
                    except ValueError:
                        # Try boolean (case-insensitive)
                        if value.lower() == 'true':
                            params[key] = True
                        elif value.lower() == 'false':
                            params[key] = False
                        else:
                            # Keep as string (strip quotes if present)
                            if (value.startswith("'") and value.endswith("'")) or \
                               (value.startswith('"') and value.endswith('"')):
                                value = value[1:-1]
                            params[key] = value

    except Exception as e:
        raise IOError(f"Error reading or parsing configuration file {filepath}: {e}")

    return params
