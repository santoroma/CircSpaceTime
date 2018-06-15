

cd //Users/gianlucamastrantonio/Dropbox/github/CircSpaceTime/src/

rm -Rf *.so
rm -Rf *.o


R CMD SHLIB   -o ProjSp.so   ProjSp.cpp -Wall-o 
R CMD SHLIB   -o ProjKrig.so   ProjKrig.cpp -Wall-o 




