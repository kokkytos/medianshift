10 DEF FNARCOS(ARG)=1.570796-ATN(ARG/SQR(1.-ARG*ARG))
20 DEF FNARCSIN(ARG)=ATN(ARG/SQR(1.-ARG*ARG))
30 DEF FNDEG(ARG)=INT(ARG)+((ARG-INT(ARG))*10.)/6.
40 DEF FNDMS(ARG)=INT(ARG)+6.*(ARG-INT(ARG))/10.
50 RD=57.29578
60 DR=1./RD
70 DIM A(4)
80 DIM B(2)
90 A(1)=-.01454
100 A(2)=-.10453
110 A(3)=-.20791
120 A(4)=.00233
130 CE=.91775
140 SE=.39715
150 INPUT "LONGITUDE IN DEG.";LO
160 INPUT "LATITUDE IN DEG.";F
170 INPUT "YEAR (4 DIGITS)";IY
180 INPUT "MONTH (NUMERAL)";IM
190 INPUT "DAY (NUMERAL)";ID
200 F=F*DR
210 C=360.
220 LI=ABS(LO)
230 SI=SIN(F)
240 CI=COS(F)
250 J=367*IY-INT(7*(IY+INT((IM+9)/12))/4)+INT(275*IM/9)+ID-730531.
260 INPUT "UNIVERSAL TIME = 0, ZONE TIME - 1, LOCAL MEAN TIME - 2";Z
270 DT=0.
280 IF Z=0. THEN LET DT=-LO/C 
290 IF Z=1. THEN LET DT=-(LI-15*INT((LI+7.5)/15))/C*SGN(LO)
300 INPUT "HOUR (4 DIGIT NUMERAL ON 24 HOUR CLOCK)";H
310 Z0=J-.5
320 IF H>0 THEN GOTO 870
330 PRINT "DATA FOR ";IY;", MONTH ";IM;", DAY ";ID
340 FOR L=1 TO 4
350 ON L GOTO 370,650,650,360
360 C=347.81
370 M=.5+DT
380 K=1
390 M=M-DT
400 E=M-LO/360.
410 GOSUB 430
420 GOTO 530
430 D=Z0+E
440 IF ABS(E)>-1 THEN LET E=E-SGN(E)
450 GOSUB 1220
460 IF L=4 THEN GOSUB 1720
470 T=T+LO+360.*E
480 T=T-INT(T/360.)*360.
490 U=T-AS
500 IF ABS(U)>180. THEN LET U=U-360.*SGN(U)
510 U=U/C
520 RETURN
530 M=M-U+DT
540 IF L<4 THEN LET K=K+1
550 ON K GOTO 600,560,600,580,600,620
560 IF M>=0. AND M<1. THEN GOTO 620
570 GOTO 590
580 IF M>=0. THEN GOTO 620
590 M=M-SGN(M)
600 K=K+1
610 GOTO 390
620 H=FNARCSIN(COS(F-DS))*RD
630 IF L-4 THEN LET H=H-.95*COS(H)
640 GOSUB 2160
650 GOSUB 2000
660 B(1)=M-H
670 B(2)=M+H
680 FOR I=1 TO 2
690 K=2*I-3
700 FOR N=1 TO 6
710 B(I)=B(I)-DT
720 E=B(I)-LO/360.
730 GOSUB 430
740 GOSUB 2000
750 B(I)=B(I)+K*H-U+DT
760 IF L<4 THEN LET N=N+1
770 ON N GOTO 820,780,820,800,820,830
780 IF B(I)>=0. AND B(I)<1. THEN GOTO 830
790 GOTO 810
800 IF B(I)>=0. THEN GOTO 830
810 B(I)=B(I)-SGN(B(I))
820 NEXT N
830 NEXT I
840 ON L GOSUB 1350,1400,1400,1610
850 NEXT L
860 END
870 INPUT "SKY CONDITION = 1,2,3,10,";SK
880 PRINT "DATA FOR ";IY;", MONTH ";IM;", DAY ";ID;", AT ";H;" HOURS"
890 E=FNDEG(H/100.)/24.-DT-LO/360.
900 D=Z0+E
910 N=1
920 GOSUB 1220
930 T=T+360.*E+LO
940 IF N=2 THEN GOSUB 1720
950 H=T-AS
960 GOSUB 2060
970 Z=H*DR
980 H=H-.95*(N-1)*COS(H*DR)
990 GOSUB 2160
1000 GOSUB 2200
1010 HA=INT(ABS(HA)+.5)*SGN(HA)
1020 ON N GOTO 1030,1090
1030 IS=133775.*M/SK
1040 PRINT "SUN AZIMUTH (DBG.) ";AZ
1050 PRINT "SUN ALTITUDE (DEG.) ";HA
1060 PRINT "SUN ILLUMINANCE (LUX) ";IS
1070 N=2
1080 GOTO 940
1090 E=FNARCOS(COS(V-LS)*CB)
1100 P=.892*EXP(-3.343/((TAN(E/2.) )^.632))+.0344*(SIN(E)-E*COS(E))
1110 P=.418*P/(1.-.005*COS(E)-.03*SIN(Z))
1120 IL=P*M/SK
1130 IS=IS+IL+.0005/SK
1140 PRINT "MOON AZIMUTH (DEG.) ";AZ
1150 PRINT "MOON ALTITUDE (DEG.) ";HA
1160 PRINT "MOON ILLUMINANCE (LUX) ";IL
1170 IL=INT(50.*(1.-COS(E))+.5)
1180 PRINT " (";IL;"% OF MOON ILLUMINATED)"
1190 PRINT "TOTAL ILLUMINANCE (LUX) ";IS
1200 GOTO 300
1210 END
1220 TD=280.46+.98565*D 
1230 T=TD-INT(TD/360)*360
1240 IF T<0. THEN LET T=T+360.
1250 TD=357.5+.9856*D
1260 G=(TD-INT(TD/360)*360)*DR
1270 LS=(T+1.91*SIN(G))*DR
1280 AS=ATN(CE*TAN(LS))*RD
1290 Y=COS(LS)
1300 IF Y<0. THEN LET AS=AS+180.
1310 SD=SE*SIN(LS)
1320 DS=FNARCSIN(SD)
1330 T=T-180.
1340 RETURN
1350 R=M
1360 GOSUB 1700
1370 PRINT "SUN MERIDIAN PASSAGE AT ";R
1380 HA=INT(ABS(HA)+.5)*SGN(HA)
1390 PRINT "ALTITUDE AT MER. PASS. ";HA
1400 FOR I=1 TO 2
1410 R=B(I)
1420 GOSUB 1700
1430 IF R>=4800. OR R<0. THEN GOTO 1680
1440 ON 2*(L-1)+I GOTO 1450,1470,1530,1550,1570,1590,1650,1670
1450 PRINT "TIME OF SUNRISE ";R
1460 GOTO 1680
1470 PRINT "TIME OF SUNSET ";R
1480 R=B(2)-B(1)
1490 IF R<0. THEN LET R=R+1.
1500 GOSUB 1700
1510 PRINT "TOTAL DAYLIGHT ";R
1520 GOTO 1680
1530 PRINT "BEGIN CIVIL TWILIGHT AT ";R
1540 GOTO 1680
1550 PRINT "END CIVIL TWILIGHT AT ";R
1550 GOTO 1680
1570 PRINT "BEGIN NAUTICAL TWILIGHT ";R
1580 GOTO 1680
1590 PRINT "END NAUTICAL TWILIGHT ";R
1600 GOTO 1680
1610 R=M
1620 GOSUB 1700
1630 PRINT "MOON MERIDIAN PASSAGE AT";R
1640 GOTO 1380
1650 PRINT "TIME OF MOONRISE ";R
1660 GOTO 1680
1670 PRINT "TIME OF MOONSET ";R
1680 NEXT I
1690 RETURN
1700 R=INT(100.*FNDMS(R*24.)+.5)
1710 RETURN
1720 TD=218.32+13.1764*D
1730 V=TD-INT(TD/360)*360
1740 IF V<0. THEN LET V=V+360.
1750 TD=134.96+13.06499*D
1760 Y=(TD-INT(TD/360)*360)*DR
1770 TD=93.27+13.22935*D
1780 O=(TD-INT(TD/360)*360)*DR
1790 TD=235.7+24.3815*D
1800 W=(TD-INT(TD/360)*360)*DR
1810 SB=SIN(Y)
1820 CB=COS(Y)
1830 X=SIN(O)
1840 S=COS(O)
1850 SD=SIN(W)
1860 CD=COS(W)
1870 V=V+(6.29-1.27*CD+.43*CB)*SB+(.66+1.27*CB)*SD-.19*SIN(G)-.23*X*S
1880 V=V*DR
1890 Y=((5.13-.17*CD)*X+(.56*SB+.17*SD)*S)*DR
1900 SV=SIN(V)
1910 SB=SIN(Y)
1920 CB=COS(Y)
1930 Q=CB*COS(V)
1940 P=CE*SV*CB-SE*SB
1950 SD=SE*SV*CB+CE*SB
1960 AS=ATN(P/Q)*RD
1970 IF Q<0. THEN LET AS=AS+180.
1980 DS=FNARCSIN(SD)
1990 RETURN
2000 H=(A(L)-SI*SD)/(CI*COS(DS))
2010 IF ABS(H)>1. THEN GOTO 2040
2020 H=FNARCOS(H)*RD/C
2030 RETURN
2040 H=1.5
2050 RETURN
2060 CD=COS(DS)
2070 CS=COS(H*DR)
2080 Q=SD*CI-CD*SI*CS
2090 P=-CD*SIN(H*DR)
2100 AZ=ATN(P/Q)*RD
2110 IF Q<0. THEN LET AZ=AZ+180.
2120 IF AZ<0. THEN LET AZ=AZ+360.
2130 AZ=INT(AZ+.5)
2140 H=FNARCSIN(SD*SI+CD*CI*CS)*RD
2150 RETURN
2160 HA=H
2170 IF H<(-5./6.) THEN GOTO 2190
2180 HA=H+1./(TAN((H+8.59/(H+4.42))*DR))/60.
2190 RETURN
2200 U=SIN(HA*DR)
2210 X=753.6616
2220 S=FNARCSIN(X*COS(HA*DR)/(X+1.))
2230 M=X*(COS(S)-U)+COS(S)
2240 M=EXP(-.21*M)*U+.0289*EXP(-.042*M)*(1.+(HA+90.)*U/57.29578)
2250 RETURN
