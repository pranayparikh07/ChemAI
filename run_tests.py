#!/usr/bin/env python
"""
ChemAI Test Runner
Convenience script to run all model tests and open results
"""

import os
import sys
import subprocess
import webbrowser
from pathlib import Path


def main():
    """Main test runner"""
    print("\n" + "=" * 80)
    print(" " * 25 + "ChemAI Model Test Runner")
    print("=" * 80)
    
    # Check if we're in the right directory
    if not os.path.exists("test_all_models.py"):
        print("\n❌ Error: test_all_models.py not found!")
        print("Please run this script from the ChemAI root directory:")
        print("  cd d:\\ChemAI")
        print("  python run_tests.py")
        return False
    
    print("\n📋 Starting comprehensive model testing...")
    print("This may take 2-5 minutes depending on your system.\n")
    
    # Run tests
    try:
        result = subprocess.run(
            [sys.executable, "test_all_models.py"],
            capture_output=False
        )
        
        if result.returncode != 0:
            print("\n⚠️  Tests completed with some warnings or errors.")
            print("Check the console output above for details.")
        else:
            print("\n✅ Tests completed successfully!")
    
    except Exception as e:
        print(f"\n❌ Error running tests: {e}")
        return False
    
    # Check if report was generated
    report_path = Path("model_test_report.html")
    if report_path.exists():
        print(f"\n✅ Report generated: {report_path.absolute()}")
        
        # Ask to open report
        print("\n" + "-" * 80)
        try:
            response = input("Would you like to open the report in your browser? (y/n): ").strip().lower()
            if response == 'y':
                webbrowser.open(str(report_path.absolute()))
                print("Opening report in browser...")
            else:
                print("You can manually open the report at:")
                print(f"  {report_path.absolute()}")
        except KeyboardInterrupt:
            print("\nSkipped opening browser.")
    else:
        print("\n⚠️  Warning: Report file not found at model_test_report.html")
    
    print("\n" + "=" * 80)
    print("Testing Complete!")
    print("=" * 80 + "\n")
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
