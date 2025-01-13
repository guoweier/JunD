# Weier Guo Python 
# Create: 12/27/2024
## usage: currently this scripts needs to be run in the folder with all fq files. 
#### UPDATE: to specify the location of fq files. #####
#### UPDATE: to write next step for threshold 1 controls=0 ####

## PACKAGES ##
import sys, math, os, time
import argparse 
import subprocess 

class RunScripts:
    
    def run_script(self, script, *args):
        command = ["python3", script]+list(args)
        print(f"Running: {' '.join(command)}")
        try:
            subprocess.run(command)
        except Exception as e:
            print(f"Failed to run {script}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Main")
    parser.add_argument("script", choices=["mapping", "binbysam", "chiread", "threshold1", "pseudojun", "threshold2", "seeds"], help="Choose the step to run.")
    parser.add_argument("args", nargs=argparse.REMAINDER, help="Arguments pass to the chosen step.")
    args = parser.parse_args()

    # map script names to their paths
    script_mapping = {
        "mapping": "../scripts/module_fqtobam.py",
        "binbysam": "../scripts/module_binbysam.py",
        "chiread": "../scripts/module_search_chiread_bin.py",
        "threshold1": "../scripts/module_threshold1_control.py",
        "pseudojun": "../scripts/module_threshold2_pseudojun.py",
        "threshold2": "../scripts/module_threshold2_sample.py",
        "seeds": "../scripts/module_collect_seeds.py"
    }

    script_path = script_mapping.get(args.script)

    # run the selected script
    runner = RunScripts()
    runner.run_script(script_path, *args.args)

if __name__ in "__main__":
    main()
