@echo off
rem Start the Reaction calculator: a local server plus the control window in your browser.
rem Extra arguments are passed on, for example:  reaction.cmd --port 8800 --no-browser
setlocal
set "ROOT=%~dp0.."
set "PYTHONPATH=%ROOT%\src;%PYTHONPATH%"
python -m reaction.app %*
