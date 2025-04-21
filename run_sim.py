# run_simulation.py

import sys, os
# Assuming parser.py is in the same directory
from src.parser import parse_input_file
from src.diatomic_gas import particle, simulation, simprops
from src.initialize import initializeGrid, initializeVelocities, initializeFiles
from src.md import md
import numpy as np

# --- Your Simulation Logic Goes Here ---
# Adapt this function to take the parameters dictionary
# and run your simulation based on its contents.
def run_simulation(params):
    """
    Runs the main simulation using parameters from the input file.

    Args:
        params (dict): A dictionary containing simulation parameters
                       read from the config file.
    """
    print("\n--- Running Simulation ---")
    print("Simulation parameters received:")
    for key, value in params.items():
        print(f"  {key}: {value} (type: {type(value).__name__})")

    # --- INTEGRATE YOUR TEST CODE HERE ---
    # Replace the placeholder logic below with your actual simulation code.
    # Access parameters from the 'params' dictionary like:
    # num_steps = params.get('num_steps', 100) # Get param or use default
    # initial_value = params.get('initial_value', 0.0)
    # model_type = params.get('model_type', 'default')
    # ... and so on for all parameters your simulation needs.

    try:
        sim = simulation()
        sim.T = params.get('T')
        sim.rho = params.get('rho')
        sim.Nm = params.get('Nm')
        sim.eq = params.get('EqSteps')
        sim.pr = params.get('PrSteps')
        sim.rc = params.get('RC')
        sim.bl = params.get('BL')
        sim.dt = params.get('dt')
        sim.outputfile = params.get('outputfile')
        sim.output = params.get('outputfreq')
        sim.k = params.get('k')
        sim.moviefile = params.get('moviefile')
        sim.moviefreq = params.get('moviefreq')
        sim.rescale_freq = params.get('rescale_freq')

        sim.itrr = 1
        sim.rdfmin = 0.8
        sim.rdfmax = 4.0
        sim.rdfN = 100
        sim.rdf = 1

        atom = initializeGrid(sim)

        initializeVelocities(sim, atom)

        initializeFiles(sim, atom)

        md(sim, atom)

    except Exception as e:
        print(f'\n Error during simulation execution: {e}')


# --- Main Execution Block ---
if __name__ == '__main__':
    # Get the configuration file path from command line arguments
    # Usage: python run_simulation.py path/to/your/simulation_config.ini
    if len(sys.argv) != 2:
        print("Usage: python run_simulation.py <config_file_path>")
        sys.exit(1)

    config_filepath = sys.argv[1]

    try:
        # Parse the configuration file
        simulation_parameters = parse_input_file(config_filepath)

        # Run the simulation with the parsed parameters
        run_simulation(simulation_parameters)

    except FileNotFoundError as e:
        print(f"Error: Configuration file not found at {config_filepath}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error parsing configuration file: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during setup or simulation: {e}", file=sys.stderr)
        sys.exit(1)
