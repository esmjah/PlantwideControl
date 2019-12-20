@echo off
title "Running Network_D1D2.genkey with C:\Program Files (x86)\Schlumberger\OLGA 2014.1.0\OlgaExecutables\OLGA-2014.1.0.exe"
pushd "C:\Git\PlantwideControl\olgaModels\NetworkOPC"
call "C:\Program Files (x86)\Schlumberger\OLGA 2014.1.0\OlgaExecutables\OLGA-2014.1.0.exe" "Network_D1D2.genkey"
title "Finished Network_D1D2.genkey"
popd
pause
exit
