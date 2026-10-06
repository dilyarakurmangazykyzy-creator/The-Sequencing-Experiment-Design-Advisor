@echo off
cd /d "%~dp0"
if not exist results\reviews mkdir results\reviews
echo Actual pipeline test run: %date% %time% > results\reviews\gaukhar_pipeline_tests.txt
python -m unittest discover -s tests -p test_pipeline.py -v >> results\reviews\gaukhar_pipeline_tests.txt 2>&1
set "gaukharCheckCode=%errorlevel%"
type results\reviews\gaukhar_pipeline_tests.txt
if not "%gaukharCheckCode%"=="0" (
  echo Checks failed. Review the saved error before reporting success.
  pause
  exit /b %gaukharCheckCode%
)
echo Tests passed. Actual output is saved in results\reviews\gaukhar_pipeline_tests.txt
echo Add your own interpretation and review notes after checking the results.
pause
