pro ta_cthph2thph,cthb,cphb,los,thb,phb

; Angles in LOS frame to vertical frame
; [rad]

thb=acos(cos(los)*cos(cthb)-cos(cphb)*sin(los)*sin(cthb))
ss=sin(cphb)*sin(cthb)/sin(thb)
cc=(cos(cthb)-cos(los)*cos(thb))/sin(los)/sin(thb)
if abs(cc) ge 1. then print,'cos(phb) is over 1: '+string(cc,format='(f12.9)') 
if ss ge 0 then phb=acos(cc>(-1.)<1.) else phb=-acos(cc>(-1.)<1.)


end
