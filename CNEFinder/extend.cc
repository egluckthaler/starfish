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
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include "cnef.h"
#include "edlib.h"

using namespace std;

bool order_qgram(QGramOcc a, QGramOcc b) 
{ 
	
	if( a.occRef == b.occRef )
	{
		return( a.occQuery < b.occQuery );
	}
	else  return ( a.occRef < b.occRef ); 

}

bool order(MimOcc a, MimOcc b) 
{ 
	
	if( a.startRef == b.startRef )
	{
		return( a.startQuery < b.startQuery );
	}
	else  return ( a.startRef < b.startRef ); 

}

bool uniqueEnt(MimOcc a, MimOcc b) 
{
	if( a.startRef == b.startRef && a.endRef == b.endRef && a.startQuery == b.startQuery && a.endQuery == b.endQuery )
	{
  		return true;
	}
	else return false;
}

int find_maximal_inexact_matches( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims, unsigned int qgram_size )
{

	sort( q_grams->begin(), q_grams->end(), order_qgram );

	fprintf ( stderr, " -Merging %i maximal exact matches\n", q_grams->size() );
	merge( sw, ref, query, q_grams, mims );

	if( mims->size() == 0 )
	{
		fprintf( stderr, "No CNEs identified!\n");
		exit(1);
	}

	q_grams->clear();

	int merged_size = sw . M * sw . l;
	fprintf ( stderr, " -Extending %i merged matches of minimum length %i, with an additional extension threshold of %.2f\n", mims->size(), merged_size, sw . s );

	#pragma omp parallel for
	for( int i=0; i<mims->size(); i++ )
	{ 	
		double minLen = min(mims->at(i).endRef-mims->at(i).startRef,mims->at(i).endQuery-mims->at(i).startQuery);
		double maxLen = max(mims->at(i).endRef-mims->at(i).startRef,mims->at(i).endQuery-mims->at(i).startQuery);

		if( mims->at(i). error / minLen < sw . t && maxLen <= sw . u )
		{	
			extend( &mims->at(i).error, (int*) &mims->at(i).startQuery, (int*) &mims->at(i).endQuery, (int*) &mims->at(i).startRef, (int*) &mims->at(i).endRef, ref, query, sw );
			adjust(  &mims->at(i).error, (int*) &mims->at(i).startQuery, (int*) &mims->at(i).endQuery, (int*) &mims->at(i).startRef, (int*) &mims->at(i).endRef, ref, query, sw );
		}
			
	}

	sort( mims->begin(), mims->end(), order );


	/* Remove overlapping CNEs */
	/*vector<MimOcc> * temp = new vector<MimOcc>;

	temp->push_back(mims->at(0));
	int i = 1;
	while( i < mims->size()  )
	{
		if(  mims->at(i).endQuery - mims->at(i).startQuery < sw . l && mims->at(i).endRef - mims->at(i).startRef < sw . l )
		{
			i++;
		}
		else if(  mims->at(i).startRef >=  temp->at( temp->size() -1 ).startRef &&  mims->at(i).endRef <=  temp->at( temp->size() -1 ).endRef &&  mims->at(i).startQuery >=  temp->at( temp->size() -1 ).startQuery &&  mims->at(i).endQuery <=  temp->at( temp->size() -1 ).endQuery )
		{
			i++;
		}
		else
		{
			if( temp->at( temp->size() -1).startRef >=  mims->at(i).startRef &&  temp->at( temp->size() -1).endRef <=  mims->at(i).endRef &&  temp->at( temp->size() -1).startQuery >=  mims->at(i).startQuery &&  temp->at( temp->size() -1).endQuery <=  mims->at(i).endQuery )	
			{
				temp->erase( temp->begin() + temp->size() -1 );
				temp->push_back( mims->at(i) );
				i++;
			}
			else
			{
				temp->push_back( mims->at(i) );
				i++;
			}
		}
	}

	mims->clear();

	for(int i=0; i<temp->size(); i++)
		mims->push_back(temp->at(i));

	delete( temp );*/


return 0;
}

int merge( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims )
{
	for( int i = 0; i<q_grams->size(); i++ )
	{	
		unsigned int current_qgram = i;	
		unsigned int edit_distance = 0;

		unsigned int q_start = q_grams->at(i).occQuery;
		unsigned int q_end = q_start + q_grams->at(i).length ;
		unsigned int r_start = q_grams->at(i).occRef;
		unsigned int r_end = r_start + q_grams->at(i).length ;
		int gap_size_ref = 0;
		int gap_size_query = 0;

		double minLen = min(r_end - r_start, q_end - q_start );
		double maxLen = max(r_end - r_start, q_end - q_start );


		for( int j = i + 1; j<q_grams->size(); j++ )
		{
		
			if( q_grams->at(j).occRef < q_grams->at(current_qgram).occRef )
				continue;

			if( maxLen >= sw . u )
				break;

			gap_size_ref = 	q_grams->at(j).occRef - ( q_grams->at(current_qgram).occRef + q_grams->at(current_qgram).length ); 
			gap_size_query = q_grams->at(j).occQuery - ( q_grams->at(current_qgram).occQuery + q_grams->at(current_qgram).length );
			

			//Check if gap in ref or query contains $ 
			bool ref$ = false;
			if( sw . p == 1 )
			{
				for(int k= q_grams->at(current_qgram).occRef + q_grams->at(current_qgram).length; k<q_grams->at(j).occRef; k++)
				{
					if( ref[k] == '$' )
					{
						ref$ = true;
						break;

					}
				}	
			}
			
			bool query$ = false;
			if( sw . p == 1 )
			{
				for(int k= q_grams->at(current_qgram).occQuery + q_grams->at(current_qgram).length ; k<q_grams->at(j).occQuery ; k++)
				{
					if( query[k] == '$' )
					{
						query$ = true;
						break;
					}
				}
			}


			if(  q_grams->at(j).occRef + q_grams->at(j).length > r_end &&  q_grams->at(j).occQuery + q_grams->at(j).length > q_end )
			{
				minLen = min(q_grams->at(j).occRef + q_grams->at(j).length - r_start, q_grams->at(j).occQuery+ q_grams->at(j).length - q_start );
				maxLen = max(q_grams->at(j).occRef + q_grams->at(j).length - r_start, q_grams->at(j).occQuery+ q_grams->at(j).length - q_start );
			}
			else 
			{	
				maxLen = maxLen;
				minLen = minLen;
			}

			if( query$ == false && ref$ == false )
			{
				if( gap_size_ref == 0 && gap_size_query / minLen <= sw . t  && gap_size_query > 0 && maxLen <= sw . u  )
				{
					if( ( edit_distance + gap_size_query )/minLen  <= sw.t  )
					{
						edit_distance = edit_distance + gap_size_query;
						q_end = q_grams->at(j).occQuery+ q_grams->at(j).length;
						r_end =  q_grams->at(j).occRef + q_grams->at(j).length;
			
						current_qgram = j;
						minLen = min(r_end - r_start, q_end - q_start );
						maxLen = max(r_end - r_start, q_end - q_start );
					}
				}
				else if( gap_size_query == 0 && gap_size_ref/minLen <= sw.t && gap_size_ref > 0 && maxLen <= sw . u) 
				{
					if( (edit_distance + gap_size_ref)/minLen <= sw.t  )
					{
						edit_distance = edit_distance + gap_size_ref;
						r_end = q_grams->at(j).occRef+ q_grams->at(j).length;
						q_end =  q_grams->at(j).occQuery + q_grams->at(j).length; 

						current_qgram = j;
						minLen = min(r_end - r_start, q_end - q_start );
						maxLen = max(r_end - r_start, q_end - q_start );
					}
				}
				else if( gap_size_query == 0 && gap_size_ref == 0  )
				{	
					r_end = q_grams->at(j).occRef + q_grams->at(j).length;
					q_end = q_grams->at(j).occQuery + q_grams->at(j).length;

					current_qgram = j;
					minLen = min(r_end - r_start, q_end - q_start );
					maxLen = max(r_end - r_start, q_end - q_start );
				}
				else if ( gap_size_query > 0 && gap_size_ref > 0 && maxLen <= sw . u)
				{	
					if( abs( gap_size_query -  gap_size_ref ) / minLen > sw.t )
						break;
				
					unsigned char * m_query = ( unsigned char * ) calloc ( gap_size_query + 1, sizeof ( unsigned char ) );
					unsigned char * m_ref = ( unsigned char * ) calloc ( gap_size_ref + 1, sizeof ( unsigned char ) );
			
					memcpy( &m_query[0], &query[ q_end ], gap_size_query );
					memcpy( &m_ref[0], &ref[ r_end ] , gap_size_ref );

					m_query[ gap_size_query ] = '\0';
					m_ref[ gap_size_ref ] = '\0';
						
					int edit_distance_temp = edit_distance + editDistanceMyers( m_query, m_ref );

					free( m_query );
					free( m_ref );

					if( edit_distance_temp/minLen <= sw.t  )
					{
						edit_distance = edit_distance_temp;
						r_end = q_grams->at(j).occRef + q_grams->at(j).length; 
						q_end = q_grams->at(j).occQuery  + q_grams->at(j).length;

						current_qgram = j;
						minLen = min(r_end - r_start, q_end - q_start );
						maxLen = max(r_end - r_start, q_end - q_start );
					}
				}	
			}
			else break;
		}



		bool longer = false;

		if( r_end-r_start >= sw . M * sw . l && r_end - r_start <= sw . u && q_end-q_start >= sw . M * sw . l && q_end - q_start <= sw . u )	
			longer = true;

		//if( r_end - r_start > q_grams->at(i).length && q_end - q_start > q_grams->at(i).length )
		//	longer = true;
		
		if ( r_end-r_start > sw . u && q_end - q_start > sw . u )
		{
			r_end = r_start + sw . l;
			q_end = q_start + sw . l;
		}
		
		if( max( r_end - r_start, q_end - q_start ) <= sw . u && longer == true )
		{
			MimOcc occ;
			occ.startRef = r_start;
			occ.endRef = r_end;
			occ.startQuery = q_start;
			occ.endQuery = q_end;
			occ.error = edit_distance;
			mims->push_back(occ);
		}

	}
	return 0;
}


int extend( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw )
{
	unsigned int toAddStartQuery = 1;
	unsigned int toAddEndQuery = 1;
	unsigned int toAddStartRef = 1;
	unsigned int toAddEndRef = 1;

	unsigned int q_start_temp = *q_start;
	unsigned int q_end_temp = *q_end;
	unsigned int r_start_temp = *r_start;
	unsigned int r_end_temp = *r_end;

	unsigned int qS = *q_start;
	unsigned int rS = *r_start;	
	unsigned int qE = *q_end;
	unsigned int rE = *r_end;

	char rsc;
	char qsc;
	char rec;
	char qec;

	unsigned int rs = *r_start;
	unsigned int qs = *q_start;
	unsigned int re = *r_end;
	unsigned int qe = *q_end;

	unsigned int edit_distance_total_L = 0;
	unsigned int edit_distance_total_R = 0;

	unsigned int edit_distance_temp = *edit_distance;
	unsigned int edit_distance_updated = *edit_distance;

	char operationEnd;
	char operationStart;

	double minLen = min( q_end_temp - q_start_temp, r_end_temp - r_start_temp );
	double maxLen = max(r_end_temp - r_start_temp, q_end_temp - q_start_temp );

	while( q_start_temp >= 0 || r_start_temp >= 0 || q_end_temp <=  strlen( ( char* ) yInput ) -1 || r_end_temp <=  strlen( ( char* ) xInput ) )
	{
		if( maxLen >= sw . u  )
			break;

		/************************************************ Score for extending right ***************************************************/
		int edit_distance_R = 0;

		char sRref;
		char sRquery;
		char iRref;
		char iRquery;
		char dRref;
		char dRquery;

		unsigned int maxSeq = max(  strlen( ( char* ) yInput ),strlen( ( char* ) xInput ) );

		if (  q_end_temp  < strlen( ( char* ) yInput )  &&  r_end_temp  < strlen( ( char* ) xInput ) ) 
		{	
			unsigned int editDist_S = 0;
			unsigned char * m_ref_R  = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_R  = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

			if( xInput[rE + toAddEndRef -1 ] == '$' || yInput[ qE + toAddEndQuery -1 ] == '$' )
			{
				editDist_S = maxSeq + 1;
			}
			else
			{
				memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
				memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
				m_ref_R[ toAddEndRef ] = '\0';
				m_query_R[ toAddEndQuery ] = '\0';

				editDist_S = editDistanceMyers( m_ref_R, m_query_R );

				sRref = m_ref_R[toAddEndRef - 1];
				sRquery = m_query_R[toAddEndQuery -1];
			}
				
			unsigned int editDist_I = 0;
			unsigned int editDist_D = 0;
			
	
			if( toAddEndRef > 1 )
			{
				if( yInput[ qE + toAddEndQuery - 1] == '$' )
					editDist_I = maxSeq + 1;
				else
				{
					memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
					memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
					m_ref_R[ toAddEndRef - 1 ] = '\0';
					m_query_R[ toAddEndQuery ] = '\0';

					editDist_I = editDistanceMyers( m_ref_R, m_query_R );

					iRref = m_ref_R[toAddEndRef - 2];
					iRquery = m_query_R[toAddEndQuery -1];
				}
	
				if( xInput[ rE + toAddEndRef - 1] == '$' )
					editDist_D = maxSeq + 1;
				else
				{
					memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef  );
					memcpy( &m_query_R[0], &yInput[qE],  toAddEndQuery  );
					m_ref_R[ toAddEndRef ] = '\0';
					m_query_R[ toAddEndQuery - 1 ] = '\0';

					editDist_D = editDistanceMyers( m_ref_R, m_query_R );

					dRref = m_ref_R[toAddEndRef - 1];
					dRquery = m_query_R[toAddEndQuery -2];
				}

			}
			else
			{
				editDist_I = maxSeq + 1;
				editDist_D = maxSeq + 1;
			}

			edit_distance_R =  min( editDist_S, min( editDist_I, editDist_D ) );

			if( edit_distance_R == editDist_S && sRref == sRquery )
			{
				operationEnd = 'S';
				rec = sRref;
				qec = sRquery;

			}
			else if( edit_distance_R == editDist_I && iRref == iRquery )
			{
				operationEnd = 'I';
				rec = iRref;
				qec = iRquery;

			}
			else if( edit_distance_R == editDist_D && dRref == dRquery )
			{
				operationEnd = 'D';
				rec = dRref;
				qec = dRquery;
			}
			else if( edit_distance_R == editDist_S )
			{
				operationEnd = 'S';
				rec = sRref;
				qec = sRquery;

			}
			else if( edit_distance_R == editDist_I  )
			{
				operationEnd = 'I';
				rec = iRref;
				qec = iRquery;
			}
			else if( edit_distance_R == editDist_D  )
			{
				operationEnd = 'D';
				rec = dRref;
				qec = dRquery;
			}

			free( m_ref_R );
			free( m_query_R );

		}
		else if( qE == strlen( ( char* ) yInput ) && rE != strlen( ( char* ) xInput ) && r_end_temp < strlen( ( char* ) xInput )  )
		{
			if( xInput[ r_end_temp + 1] == '$' )
				edit_distance_R = maxSeq + 1;
			else
			{
				edit_distance_R = edit_distance_total_R + 1;
				operationEnd = 'D';

				rec = xInput[ r_end_temp  ];
				qec = yInput[ strlen( ( char* ) yInput ) - 1];
			}
		}
		else if( rE == strlen( ( char* ) xInput ) && qE != strlen( ( char* ) yInput ) && q_end_temp < strlen( ( char* ) yInput ) )
		{
			
			if( yInput[ q_end_temp + 1] == '$' )
				edit_distance_R = maxSeq + 1;
			else
			{
				edit_distance_R = edit_distance_total_R + 1;
				operationEnd = 'I';

				qec = yInput[ q_end_temp  ];
				rec = xInput[ strlen( ( char* ) xInput ) - 1];
			}
		}
		else if ( q_end_temp  < strlen( ( char* ) yInput ) && r_end_temp >= strlen( ( char* ) xInput ) )	
		{
		
			if( yInput[ qE + toAddEndQuery  -1 ] == '$' )
				edit_distance_R = maxSeq + 1;
			else
			{
				unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
				unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

				memcpy( &m_ref_R[0], &xInput[rE],  strlen( ( char* ) xInput ) - rE );
				memcpy( &m_query_R[0], &yInput[qE], toAddEndQuery  );
				m_ref_R[ toAddEndRef ] = '\0';
				m_query_R[ toAddEndQuery ] = '\0';
		
				edit_distance_R =  editDistanceMyers( m_ref_R, m_query_R );
				operationEnd = 'I';

				rec = m_ref_R[ strlen( ( char* ) m_ref_R )  - 1];
				qec = m_query_R[ toAddEndQuery - 1];

				free( m_ref_R );
				free( m_query_R );
			}
		}
		else if ( q_end_temp  >= strlen( ( char* ) yInput ) - 1 && r_end_temp < strlen( ( char* ) xInput ) - 1 )	
		{
			if( xInput[ rE+toAddEndRef  -1 ] == '$' )
				edit_distance_R = maxSeq + 1;
			else
			{
				unsigned char * m_ref_R = ( unsigned char * ) calloc (  toAddEndRef + 1, sizeof ( unsigned char ) );
				unsigned char * m_query_R = ( unsigned char * ) calloc ( toAddEndQuery + 1, sizeof ( unsigned char ) );

				memcpy( &m_ref_R[0], &xInput[rE],  toAddEndRef );
				memcpy( &m_query_R[0], &yInput[qE], strlen( ( char* ) yInput ) - qE   );
				m_ref_R[ toAddEndRef ] = '\0';
				m_query_R[ toAddEndQuery ] = '\0';
		
				edit_distance_R =  editDistanceMyers( m_ref_R, m_query_R );
				operationEnd = 'D';

				rec = m_ref_R[ toAddEndRef - 1 ];
				qec = m_query_R[  strlen( ( char* ) m_query_R )  - 1];

				free( m_ref_R );
				free( m_query_R );
			}
		}
		else 
		{	
			edit_distance_R = maxSeq + 1;
			rec = xInput[ strlen( ( char * ) xInput ) - 1 ];
			qec = yInput[ strlen( ( char * ) yInput ) - 1 ];

		}


		/*********************************************** score for extending left *************************************************/
		unsigned int edit_distance_L = 0;

		char sLref;
		char sLquery;
		char iLref;
		char iLquery;
		char dLref;
		char dLquery;

		if(  q_start_temp  > 0 &&  r_start_temp > 0   )  
		{
			unsigned int editDist_S;
			unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
			unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
			

			if( xInput [rS - toAddStartRef] == '$' || yInput [qS - toAddStartQuery] == '$' )
			{
				editDist_S  = maxSeq +1;
			}	
			else
			{

				memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
				memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
				m_ref_L[ toAddStartRef ] = '\0';
				m_query_L[ toAddStartQuery ] = '\0';

				editDist_S = editDistanceMyers( m_ref_L, m_query_L );

				sLref = m_ref_L[0];
				sLquery = m_query_L[0];
			}

			unsigned int editDist_I = 0;
			unsigned int editDist_D = 0;

			if( toAddStartRef > 1 )
			{
				if(  yInput [qS - toAddStartQuery] == '$' )
				{
					editDist_I  = maxSeq +1;
				}
				else
				{
					memcpy( &m_ref_L[0], &xInput [ rS - toAddStartRef + 1  ], toAddStartRef  );
					memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
					m_ref_L[ toAddStartRef - 1 ] = '\0';
					m_query_L[ toAddStartQuery ] = '\0';
					
					editDist_I = editDistanceMyers( m_ref_L, m_query_L );
				
					iLref = m_ref_L[0];
					iLquery = m_query_L[0];
				}


				if( xInput [rS - toAddStartRef] == '$' )
					editDist_D  = maxSeq +1;
				else
				{
					memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
					memcpy( &m_query_L[0], &yInput [ qS - toAddStartQuery + 1 ], toAddStartQuery   );
					m_ref_L[ toAddStartRef ] = '\0';
					m_query_L[ toAddStartQuery - 1 ] = '\0';

					editDist_D = editDistanceMyers( m_ref_L, m_query_L );
	
					dLref = m_ref_L[0];
					dLquery = m_query_L[0];
				}
			}
			else
			{
				editDist_I = maxSeq + 1;
				editDist_D = maxSeq + 1;
			}

			
			edit_distance_L =  min( editDist_S, min( editDist_I, editDist_D ) );

			if( edit_distance_L == editDist_S && sLref == sLquery  )
			{
				operationStart = 'S';
				rsc = sLref;
				qsc = sLquery;
			}
			else if( edit_distance_L == editDist_I && iLref == iLquery   )
			{
				operationStart = 'I';
				rsc = iLref;
				qsc = iLquery;
			}
			else if( edit_distance_L == editDist_D && dLref == dLquery  )
			{
				operationStart = 'D';
				rsc = dLref;
				qsc = dLquery;
			}
			else if( edit_distance_L == editDist_S  )
			{
				operationStart = 'S';
				rsc = sLref;
				qsc = sLquery;
			}
			else if( edit_distance_L == editDist_I  )
			{
				operationStart = 'I';
				rsc = iLref;
				qsc = iLquery;
			}
			else if( edit_distance_L == editDist_D  )
			{
				operationStart = 'D';
				rsc = dLref;
				qsc = dLquery;
			}

			free( m_ref_L );
			free( m_query_L );

		}
		else if( qS == 0 && rS != 0 && r_start_temp > 0 )
		{
			if( xInput [r_start_temp - 1] == '$'  )
			{	
				edit_distance_L  = maxSeq +1;
			}
			else
			{
				edit_distance_L = edit_distance_total_L + 1;
				operationStart = 'D';

				rsc = xInput[ r_start_temp - 1 ];
				qsc = yInput[0];
			}
		}
		else if( rS == 0 && qS != 0 && q_start_temp > 0 )
		{	
			if( yInput [q_start_temp - 1] == '$'  )
			{	
				edit_distance_L  = maxSeq +1;
			}
			else
			{
				edit_distance_L = edit_distance_total_L + 1;
				operationStart = 'I';

				qsc = yInput[ q_start_temp - 1];
				rsc = xInput[0];
			}
		}
		else if ( q_start_temp  <= 0 && r_start_temp > 0 )	
		{
			if( xInput [rS - toAddStartRef] == '$'  )
			{
				edit_distance_L  = maxSeq +1;
			}
			else
			{
				unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
				unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
				
				memcpy( &m_ref_L[0], &xInput [rS - toAddStartRef], toAddStartRef );
				memcpy( &m_query_L[0], &yInput [0], qS );
				m_ref_L[ toAddStartRef ] = '\0';
				m_query_L[ toAddStartQuery ] = '\0';
			
				edit_distance_L = editDistanceMyers( m_ref_L, m_query_L );
				operationStart = 'D';

				rsc = m_ref_L[0];
				qsc = m_query_L[0];

				free( m_ref_L );
				free( m_query_L );
			}

		}
		else if ( q_start_temp  > 0 && r_start_temp <= 0 )	
		{	
			if( yInput [qS - toAddStartQuery] == '$' )
				edit_distance_L = maxSeq + 1;
			else
			{
				unsigned char * m_ref_L = ( unsigned char * ) calloc ( toAddStartRef + 1, sizeof ( unsigned char ) );
				unsigned char * m_query_L = ( unsigned char * ) calloc ( toAddStartQuery + 1, sizeof ( unsigned char ) );
				
				memcpy( &m_ref_L[0], &xInput [0], rS );
				memcpy( &m_query_L[0], &yInput [qS - toAddStartQuery], toAddStartQuery );
				m_ref_L[ toAddStartRef ] = '\0';
				m_query_L[ toAddStartQuery ] = '\0';

				edit_distance_L = editDistanceMyers( m_ref_L, m_query_L );
				operationStart = 'I';

				rsc = m_ref_L[0];
				qsc = m_query_L[0];

				free( m_ref_L );
				free( m_query_L );

			}
		}
		else 
		{	
			rsc = xInput[0];
			qsc = yInput[0];
			edit_distance_L = maxSeq + 1;
		}

		
		/*********************************************** computing extension *************************************************/
		int direction = 1;
		if(  (edit_distance_L + edit_distance_R + edit_distance_temp) / (minLen+1) > sw .t )
		{
			if(  edit_distance_temp + edit_distance_total_R + edit_distance_L < edit_distance_temp + edit_distance_R + edit_distance_total_L  && (edit_distance_temp + edit_distance_total_R + edit_distance_L ) / ( minLen+1 ) <= sw.t + sw.s ) //extend left
			{
				if( operationStart == 'S' )
				{
					q_start_temp--;
					r_start_temp--;
					toAddStartQuery++;
					toAddStartRef++;
					
				}
				else if( operationStart == 'I' )
				{
					toAddStartQuery++;
					q_start_temp--;
				}

				else if( operationStart == 'D' )
				{
					toAddStartRef++;
					r_start_temp--;
				}

			
				edit_distance_total_L = edit_distance_L;
				edit_distance_updated = edit_distance_temp + edit_distance_total_R + edit_distance_L;

				minLen = min( q_end_temp - q_start_temp, r_end_temp - r_start_temp );
				maxLen = max( q_end_temp - q_start_temp, r_end_temp - r_start_temp );

				if( rsc == qsc && ( edit_distance_updated / minLen ) <= sw.t  )
				{
			
					rs = r_start_temp;
					qs = q_start_temp;
				}
				direction = 0;
				
		
			}
			else if (  edit_distance_temp + edit_distance_R + edit_distance_total_L < edit_distance_temp + edit_distance_total_R + edit_distance_L && (edit_distance_temp + edit_distance_R + edit_distance_total_L ) / ( minLen+1 ) <= sw .t + sw.s) //extend right
			{
				if( operationEnd == 'S' )
				{
					q_end_temp++;
					r_end_temp++;
					toAddEndQuery++;
					toAddEndRef++;
				}
				else if( operationEnd == 'I' )
				{		
					toAddEndQuery++;

					q_end_temp++;
				}
				else if( operationEnd == 'D' )
				{	
					toAddEndRef++;
					r_end_temp++;
				}

				edit_distance_total_R = edit_distance_R;
				edit_distance_updated = edit_distance_temp + edit_distance_R + edit_distance_total_L;

				minLen = min( q_end_temp - q_start_temp, r_end_temp - r_start_temp );
				maxLen = max( q_end_temp - q_start_temp, r_end_temp - r_start_temp );

				if( rec == qec && edit_distance_updated/minLen <= sw.t)
				{
					re = r_end_temp;
					qe = q_end_temp;
				}
				direction = 1;
				
			}

			else if (  edit_distance_temp + edit_distance_R + edit_distance_total_L == edit_distance_temp + edit_distance_total_R + edit_distance_L && ( edit_distance_temp + edit_distance_R + edit_distance_total_L )/ ( minLen + 1 ) <=sw.t + sw.s) //extend based on previous extension
			{
				if( direction == 0 )
				{
					if( operationStart == 'S' )
					{
						q_start_temp--;
						r_start_temp--;
						toAddStartQuery++;
						toAddStartRef++;
					
					}
					else if( operationStart == 'I' )
					{
						toAddStartQuery++;
						q_start_temp--;
					}

					else if( operationStart == 'D' )
					{
						toAddStartRef++;
						r_start_temp--;
					}

			
					edit_distance_total_L = edit_distance_L;
					edit_distance_updated = edit_distance_temp + edit_distance_total_R + edit_distance_L;

					minLen = min( q_end_temp - q_start_temp, r_end_temp - r_start_temp );	
					maxLen = max( q_end_temp - q_start_temp, r_end_temp - r_start_temp );

					if( rsc == qsc && edit_distance_updated/minLen <=sw.t )
					{
			
						rs = r_start_temp;
						qs = q_start_temp;
					}
					direction = 0;
				}
				else if( direction == 1 )
				{	
					if( operationEnd == 'S' )
					{
						q_end_temp++;
						r_end_temp++;
						toAddEndQuery++;
						toAddEndRef++;
					}
					else if( operationEnd == 'I' )
					{		
						toAddEndQuery++;

						q_end_temp++;
					}
					else if( operationEnd == 'D' )
					{	
						toAddEndRef++;
						r_end_temp++;
					}

					edit_distance_total_R = edit_distance_R;
					edit_distance_updated = edit_distance_temp + edit_distance_R + edit_distance_total_L;

					minLen = min( q_end_temp - q_start_temp, r_end_temp - r_start_temp );
					maxLen = max( q_end_temp - q_start_temp, r_end_temp - r_start_temp );

					if( rec == qec && edit_distance_updated/minLen <= sw.t)
					{
						re = r_end_temp;
						qe = q_end_temp;
					}
					direction = 1;	
				}
				
			}
			else break;
			
		}
		else if( ( edit_distance_temp +  edit_distance_L + edit_distance_R)  / (minLen + 2) <= sw.t ) //extend both directions
		{ 

			if( operationEnd == 'S' )
			{
				q_end_temp++;
				r_end_temp++;
				toAddEndQuery++;
				toAddEndRef++;
			}
			else if( operationEnd == 'I' )
			{	
				toAddEndQuery++;
				q_end_temp++;
			}

			else if( operationEnd == 'D' )
			{
				toAddEndRef++;
				r_end_temp++;
			}
			
			if( operationStart == 'S' )
			{
				q_start_temp--;
				r_start_temp--;
				toAddStartRef++;
				toAddStartQuery++;

			}
			else if( operationStart == 'I' )
			{
				
				toAddStartQuery++;
				q_start_temp--;
			}
			
			else if( operationStart == 'D' )
			{
				r_start_temp--;
				toAddStartRef++;

			}

			edit_distance_total_L = edit_distance_L;
			edit_distance_total_R = edit_distance_R;

			edit_distance_updated =  edit_distance_temp + edit_distance_L + edit_distance_R;
			minLen = min( q_end_temp - q_start_temp, r_end_temp - r_start_temp );
			maxLen = max( q_end_temp - q_start_temp, r_end_temp - r_start_temp );

			if( rsc == qsc && edit_distance_updated/minLen <= sw.t)
			{
			
				rs = r_start_temp;
				qs = q_start_temp;
			}
		
			if( rec == qec && edit_distance_updated/minLen <= sw.t )
			{
				re = r_end_temp;
				qe = q_end_temp;
			}

			
		}
	}


	if( qs < 0 )
		*q_start = 0;
	else *q_start = qs;
		
	if( rs < 0 )
		*r_start = 0;
	else *r_start = rs;

	*edit_distance = edit_distance_updated;

	
	if( qe > strlen( ( char* ) yInput ) )
		*q_end =  strlen( ( char* ) yInput );
	else *q_end = qe;

	if( re >  strlen( ( char* ) xInput ) )
		*r_end = strlen( ( char* ) xInput );
	else *r_end = re;
	
return 0;
}

int adjust( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw )
{
	unsigned int rS = *r_start;
	unsigned int qS = *q_start;
	unsigned int rE = *r_end;
	unsigned int qE = *q_end;

	unsigned int rSb = *r_start;
	unsigned int rEb = *r_end;
	unsigned int qSb = *q_start;
	unsigned int qEb = *q_end;

	unsigned char * A = ( unsigned char * ) calloc (  rE - rS + 1, sizeof ( unsigned char ) );
	unsigned char * B = ( unsigned char * ) calloc ( qE - qS + 1, sizeof ( unsigned char ) );

	memcpy( &A[0], &xInput[rS],  rE - rS  );
	memcpy( &B[0], &yInput[qS],  qE - qS);
	A[ rE - rS ] = '\0';
	B[ qE- qS ] = '\0';
		
        *edit_distance = editDistanceMyers( A, B );

	free( A );
	free( B );

	
	unsigned int minLen = min( rE - rS, qE - qS );
	unsigned int maxLen = max( rE - rS, qE - qS );
	
	if(  *edit_distance / minLen < sw.t && maxLen <= sw . u  ) 
	{
		unsigned int eD = *edit_distance;

		extend( ( unsigned int*) &eD, (int*) &qS, (int*) &qE, (int*) &rS, (int*) &rE,  xInput, yInput, sw );

		*q_start = qS;
		*q_end = qE;
		*r_start = rS;
		*r_end = rE;

		unsigned char * A2 = ( unsigned char * ) calloc (  rE - rS + 1, sizeof ( unsigned char ) );
		unsigned char * B2 = ( unsigned char * ) calloc ( qE - qS + 1, sizeof ( unsigned char ) );

		memcpy( &A2[0], &xInput[rS],  rE - rS  );
		memcpy( &B2[0], &yInput[qS],  qE - qS );
		A2[ rE - rS ] = '\0';
		B2[ qE- qS ] = '\0';
		
		*edit_distance =  editDistanceMyers( A2, B2 );

		free( A2 );
		free( B2 );

	}

	while( rSb != *r_start || rEb != *r_end || qSb != *q_start || qEb != *q_end )
	{
		unsigned int eD = *edit_distance;
		
		rSb = *r_start;
		rEb = *r_end;
		qSb = *q_start;
		qEb = *q_end;

		extend( ( unsigned int*) &eD, (int*) &qS, (int*) &qE, (int*) &rS, (int*) &rE,  xInput, yInput, sw );

		*q_start = qS;
		*q_end = qE;
		*r_start = rS;
		*r_end = rE;

		unsigned char * A2 = ( unsigned char * ) calloc (  rE - rS + 1, sizeof ( unsigned char ) );
		unsigned char * B2 = ( unsigned char * ) calloc ( qE - qS + 1, sizeof ( unsigned char ) );

		memcpy( &A2[0], &xInput[rS],  rE - rS  );
		memcpy( &B2[0], &yInput[qS],  qE - qS );
		A2[ rE - rS ] = '\0';
		B2[ qE- qS ] = '\0';
		
		*edit_distance =  editDistanceMyers( A2, B2 );

		free( A2 );
		free( B2 );
	}
	
return 0;
}

/*
Myers Bit-Vector algorithm implemented using edlib Library
*/
int editDistanceMyers( unsigned char * xInput, unsigned char * yInput )
{
	unsigned int score = edlibAlign( (const char*) xInput, strlen( (char*) xInput ), (const char*) yInput, strlen( (char*) yInput ), edlibDefaultAlignConfig()).editDistance;

	return score;
}
