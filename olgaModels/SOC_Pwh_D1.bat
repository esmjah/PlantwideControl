@echo off
title "Running SOC_Pwh_D1.genkey with C:\Program Files (x86)\Schlumberger\OLGA 2014.3.0\OlgaExecutables\OLGA-2014.3.0.exe"
pushd "C:\olgaSOC2014"
call "C:\Program Files (x86)\Schlumberger\OLGA 2014.3.0\OlgaExecutables\OLGA-2014.3.0.exe" "SOC_Pwh_D1.genkey"
title "Finished SOC_Pwh_D1.genkey"
popd
pause
exit
