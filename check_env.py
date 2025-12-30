import sys
import importlib.util

REQUIRED_PYTHON = (3, 8)

REQUIRED_PACKAGES = [
    "numpy",
    "scipy",
    "matplotlib",
    "pandas",
    "sklearn",
    "cv2",
    "PIL",
    "skimage",
    "shapely",
    "IPython",
    "cigsegy"
    #"TkinterDnD2",
]

def check_python():
    if sys.version_info < REQUIRED_PYTHON:
        print(
            f"Python {REQUIRED_PYTHON[0]}.{REQUIRED_PYTHON[1]}+ required.\n"
            f"Detected: {sys.version.split()[0]}"
        )
        return False
    return True

def check_packages():
    missing = []
    for pkg in REQUIRED_PACKAGES:
        if importlib.util.find_spec(pkg) is None:
            missing.append(pkg)
    return missing


def check_environment():
    if not check_python():
        return ["<python_version>"]

    return check_packages()

if __name__ == "__main__":
    print("Checking environment...\n")

    missing = check_environment()

    if missing:
        print("Missing requirements:")
        for m in missing:
            print("  -", m)

        print("\nPlease run:")
        print("  pip install -r requirements.txt")

        print("\nIf problems persist, see:")
        print("  environment_info.txt")
        sys.exit(1)

    print("Environment OK ✅")
