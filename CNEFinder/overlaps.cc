/**
    CNEFinder
    Copyright (C) 2017 Lorraine A. K. Ayad, Solon P. Pissis, Dimitris Polychronopoulos

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <sys/time.h>
#include <algorithm>
#include <omp.h>
#include "cnef.h"

using namespace std;


int remove_overlaps( vector<MimOcc> * mims, TSwitch sw )
{
	vector<MimOcc> * temp = new vector<MimOcc>;

	temp->push_back(mims->at(0));
	int i = 1;
	while( i < mims->size()  )
	{

		
		int mimEndQuery = mims->at(i).endQuery;
		int mimStartQuery = mims->at(i).startQuery;
		int mimEndRef = mims->at(i).endRef;
		int mimStartRef = mims->at(i).startRef;
		
		int tempEndQuery =  temp->at( temp->size() -1 ).endQuery;
		int tempStartQuery = temp->at( temp->size() -1).startQuery;
		int tempEndRef =  temp->at( temp->size() -1 ).endRef;
		int tempStartRef =  temp->at( temp->size() -1).startRef;

		if(  mimEndQuery - mimStartQuery < sw . l && mimEndRef - mimStartRef < sw . l )
		{
			i++;
		}
		else if(  mimStartRef >=  tempStartRef &&  mimEndRef <=  tempEndRef &&  mimStartQuery >= tempStartQuery &&  mimEndQuery <=  tempEndQuery )
		{ 
			i++;
		}
		else if( tempStartRef >= mimStartRef && tempEndRef <=  mimEndRef &&  tempStartQuery >=  mimStartQuery &&  tempEndQuery <= mimEndQuery )	
		{
			temp->erase( temp->begin() + temp->size() -1 );
			temp->push_back( mims->at(i) );
			i++;
		}	
		else if( tempStartRef <= mimStartRef && tempEndRef >=  mimEndRef &&  mimStartQuery >= tempStartQuery && mimStartQuery < tempEndQuery && mimEndQuery >= tempEndQuery  || tempStartRef <= mimStartRef && tempEndRef >=  mimEndRef  && mimStartQuery <= tempStartQuery && tempStartQuery < mimEndQuery && mimEndQuery <= tempEndQuery )
		{
			i++;

		}
		else if ( tempStartQuery <=  mimStartQuery &&  tempEndQuery >= mimEndQuery &&  mimStartRef >= tempStartRef && mimStartRef < tempEndRef && mimEndRef >= tempEndRef  || tempStartQuery <=  mimStartQuery &&  tempEndQuery >= mimEndQuery && mimStartRef >= tempStartRef && mimStartRef < tempEndRef && mimEndRef >= tempEndRef )
		{
			i++;
		}
		else if( mimStartRef >= tempStartRef && mimStartRef < tempEndRef && mimEndRef >= tempEndRef && mimStartQuery >= tempStartQuery && mimStartQuery < tempEndQuery && mimEndQuery >= tempEndQuery || mimStartRef >= tempStartRef && mimStartRef < tempEndRef && mimEndRef >= tempEndRef && mimStartQuery <= tempStartQuery && tempStartQuery < mimEndQuery && mimEndQuery <= tempEndQuery  )	
		{	
			int lenMimRef = mimEndRef - mimStartRef;
			int lenTempRef = tempEndRef - tempStartRef;

			int lenMimQuery = mimEndQuery - mimStartQuery;
			int lenTempQuery = tempEndQuery  - mimStartQuery;

				
			if( lenMimRef > lenTempRef && lenMimQuery > lenTempRef )
			{
				temp->erase( temp->begin() + temp->size() -1 );
				temp->push_back( mims->at(i) );
				i++;
			}
			else if( lenMimRef >= lenTempRef && lenMimQuery <= lenTempQuery )
			{
				if(  lenMimRef - lenTempRef > lenTempQuery - lenMimQuery )
				{	
					temp->erase( temp->begin() + temp->size() -1 );
					temp->push_back( mims->at(i) );
						
				}
				i++;
			}
			else if( lenMimRef <= lenTempRef && lenMimQuery <= lenTempRef )
			{
				i++;
			}
			else if( lenMimRef <= lenTempRef && lenMimQuery >= lenTempQuery )
			{
				if(  lenTempRef-lenMimRef < lenMimQuery - lenTempQuery )
				{	
					temp->erase( temp->begin() + temp->size() -1 );
					temp->push_back( mims->at(i) );
						
				}
				i++;
			}
			else i++;

		}
		else
		{
			temp->push_back( mims->at(i) );
			i++;
		}
	}

	mims->clear();

	for(int i=0; i<temp->size(); i++)
		mims->push_back(temp->at(i));

	delete( temp );

return 0;
}

