@echo off
title "Running Network_Nominal.genkey with C:\Program Files (x86)\Schlumberger\OLGA 2014.1.0\OlgaExecutables\OLGA-2014.1.0.exe"
pushd "C:\Git\PlantwideControl\olgaModels\NetworkOPC"
call "C:\Program Files (x86)\Schlumberger\OLGA 2014.1.0\OlgaExecutables\OLGA-2014.1.0.exe" "Network_Nominal.genkey"
title "Finished Network_Nominal.genkey"
popd
pause
exit
