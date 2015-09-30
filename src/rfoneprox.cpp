//----------------------------------------------------------------------
// File:                rfoneprox.cpp
// Programmer:          Nicholas Crookston
// Description:         computes a proximity vector for one observation
// Last modified:       2008/11/18 (Version 1.0: uses a binary search approach)
//----------------------------------------------------------------------
// Copyright Notice:
// This code was written and prepared by a U.S. Government employee on official 
// time and therefore is in the public domain and not subject to copyright.
//
//----------------------------------------------------------------------
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

extern "C" {
 SEXP rfoneprox(SEXP nodes, SEXP srt, SEXP nobs, SEXP ntree, SEXP obs, SEXP prox)
 {
   int itr,myNode,mid,top,bot,me,it,ptr;   
   for (itr=0; itr<INTEGER(ntree)[0]; itr++) 
   {
     myNode=INTEGER(obs)[itr];
     top=0;
     bot=INTEGER(nobs)[0];
     it=-1;
     while (it==-1)
     {
       if (bot-top < 6) 
       {
         for (mid=top; mid<bot; mid++) 
         {
           ptr = INTEGER(srt)[INTEGER(nobs)[0]*itr+mid];
           me  = INTEGER(nodes)[INTEGER(nobs)[0]*itr+ptr];
           if (myNode == me) 
           {
             it=mid; 
             break;
           }
         }
         // the node has no match in the training data for the current tree.
         if (it==-1) break;  
       } 
       else
       {   
         mid=(top+bot)/2;
         ptr = INTEGER(srt)[INTEGER(nobs)[0]*itr+mid];
         me  = INTEGER(nodes)[INTEGER(nobs)[0]*itr+ptr];
         if (myNode == me) 
         {
           it=mid;
           break;
         }
         else if (myNode > me) top=mid;
         else bot=mid;
       }
     }   
     if (it>-1) 
     {     
       for (mid=it; mid>=0; mid--)
       {
         ptr = INTEGER(srt)[INTEGER(nobs)[0]*itr+mid];
         me  = INTEGER(nodes)[INTEGER(nobs)[0]*itr+ptr]; 
         if (myNode != me) break;
         INTEGER(prox)[ptr]++;
       }
       for (mid=it+1; mid<bot; mid++)
       {
         ptr = INTEGER(srt)[INTEGER(nobs)[0]*itr+mid];
         me  = INTEGER(nodes)[INTEGER(nobs)[0]*itr+ptr];
         if (myNode != me) break;
         INTEGER(prox)[ptr]++;
       }
     }
   }
   return(prox);
 }
}
