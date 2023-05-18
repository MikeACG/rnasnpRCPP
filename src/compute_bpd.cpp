#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_bpd(NumericMatrix wild, NumericMatrix mut, int regionX, int regionY) {

    int startCompute = 0;
    int endCompute = wild.nrow() - 1;
    List result = List::create();

    // Function to compute distance using RNAplfold matrix 
    int pos, posI,posJ;
    int localityX, localityY;
    localityX = regionX; localityY=regionY;
    double sum;double max;int maxpos;
    double previousSum;
    max=0;
    sum=0;

    //Compute first element for the recursion computation of dk
    pos=startCompute;
    for(posI=pos; posI < pos+localityX && posI <= endCompute; posI++){
        for(posJ=posI; posJ < posI + localityY && posJ <= endCompute; posJ++){
            double diff;
            diff=wild(posI, posJ) - mut(posI, posJ);
            sum+=diff*diff;
        }
    }

    //##############################
    //Start recursion for 5UTR
    //##############################
    maxpos=startCompute;
    max=sum;
    for(pos=startCompute+1; pos<= endCompute; pos++){
        posI = pos + localityX+1;
        for(posJ=posI; posJ < posI + localityY && posJ <= endCompute; posJ++){
            double diff;
            diff=wild(posI, posJ) - mut(posI, posJ);
            sum+=diff*diff;
        }
        posI = pos;
        for(posJ=posI; posJ < posI + localityY && posJ <= endCompute; posJ++){
            double diff;
            diff=wild(posI, posJ) - mut(posI, posJ);
            sum-=diff*diff;
        }
        if(max<sum){
            max=sum;
            maxpos=pos;
        }
    }
    result.push_back(max, "dmax");
    result.push_back(maxpos, "start");

    //#################################
    // Function to find the exact interval of max difference
    //##################################

    int maxend;max=0;
    for(posJ=maxpos+localityX;posJ<maxpos+localityX+localityY && posJ<=endCompute;posJ++){
        sum=0;
        
        for(posI=maxpos;posI<posJ;posI++){
            for(pos=(posI+1);pos<=posJ;pos++)
            {
                double diff;
                diff = wild(posI, pos) - mut(posI, pos);
                sum+=diff*diff;
            }
        }
        //printf("%d %d\n",maxpos,posJ);
        sum/=(posJ-maxpos);
        if(max<sum){
            max=sum;
            maxend=posJ;
        }
    }
    result.push_back(max, "d");
    result.push_back(maxend, "end");

    return result;
}
