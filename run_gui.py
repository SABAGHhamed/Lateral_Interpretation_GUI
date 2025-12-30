import sys
import ctypes

mutex = ctypes.windll.kernel32.CreateMutexW(None, False, "LateralInterpretationGUI_MUTEX")

from check_env import check_packages
from seismic_gui.app import main

if ctypes.windll.kernel32.GetLastError() == 183:
    # ERROR_ALREADY_EXISTS
    sys.exit(0)


missing = check_packages()
if missing:
    raise RuntimeError(
        "Missing dependencies. Run:\n"
        "  pip install -r requirements.txt"
    )

main()
