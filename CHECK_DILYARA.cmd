@echo off
cd /d "%~dp0"
if not exist results\reviews mkdir results\reviews
echo Actual advisor test run: %date% %time% > results\reviews\dilyara_advisor_tests.txt
python -m unittest discover -s tests -p test_advisor.py -v >> results\reviews\dilyara_advisor_tests.txt 2>&1
set "dilyaraCheckCode=%errorlevel%"
type results\reviews\dilyara_advisor_tests.txt
if not "%dilyaraCheckCode%"=="0" (
  echo Checks failed. Review the saved error before reporting success.
  pause
  exit /b %dilyaraCheckCode%
)
echo Tests passed. Actual output is saved in results\reviews\dilyara_advisor_tests.txt
echo Add your own interpretation and scenario notes after checking the results.
pause
