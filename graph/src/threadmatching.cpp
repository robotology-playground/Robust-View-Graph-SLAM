/*
 *  Created on: Nov 11, 2016
 *      Author: Nicolo' Genesio
 *      email: nicolo.genesio@iit.it
 */

#include "threadmatching.h"



ThreadMatching::ThreadMatching(vgSLAMBuffer<SlamType> &bufferIn, vgSLAMBuffer<SlamType> &bufferOut,
                               cv::Ptr<cv::DescriptorMatcher> _matcher): vgSLAMThread(bufferIn, bufferOut),matcher(_matcher){
    first = second = NULL;

}
ThreadMatching::~ThreadMatching(){
    if(first) {
        first->free();
        delete first;
    }

    if(second) {
        second->free();
        delete second;
    }
}

void ThreadMatching::run (){

    while(!interrupted) {
        yInfo()<<"ThreadMatching:Reading buffer";
        SlamType* data = new SlamType();
        MatchesVector mv;
        if(!bufferIn->read(*data)) {
            delete data;
            continue;
        }

        if(!first) {
            first = data;
            continue;
        }

        if(!second) {
            second = data;
            data->matching=new MatchesVector;
            matcher->match(*first->descriptor,*second->descriptor,mv,cv::noArray());
            yDebug()<<"ThreadMatching first->second: size matching"<<mv.size();
            continue;
        }

        // processing
        yDebug()<<"ThreadMatching: first:"<<first->stamp->getCount()<<"second:"<<second->stamp->getCount()<<"third:"<<data->stamp->getCount();//tested ok
        // ...
        matcher->match(*second->descriptor,*data->descriptor,mv,cv::noArray());
        yDebug()<<"ThreadMatching second->third: size matching"<<mv.size();
        matcher->match(*first->descriptor,*data->descriptor,mv,cv::noArray());
        yDebug()<<"ThreadMatching first->third: size matching"<<mv.size();


        //substitution
        first->free();
        delete first;
        first = second;
        second = data;

        yDebug()<<"ThreadMatching:count="<<countProcessed;

        countProcessed++;

    } //end while
}
