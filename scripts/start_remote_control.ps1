# start_remote_control.ps1 -- IKMHIWI03: start `claude remote-control` in its own PowerShell window, so it keeps
# running when the VS Code tunnel terminal or the phone's browser is closed (21 Oct 2026). The session then shows up
# in the Claude Code app; give it instructions there. ANSYS chains run detached (run_chain.ps1) either way.
#   .\scripts\start_remote_control.ps1
param([string]$Repo = "C:\Users\nishioka\git\pde-fem-biofilm")
Start-Process powershell -ArgumentList '-NoExit', '-Command', "cd $Repo; claude remote-control"
