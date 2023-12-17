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

#include <vector>
#include <string>
#define INITIAL_SC		-100000
#define ALLOC_SIZE               104857
#define NA			'N'
#define GAP 			'-'
#define DOL			'$'
#define INS			1
#define DEL			1
#define SUB			1
#define MAT			0

#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)  
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))

using namespace std;

struct TSwitch
 {
   char               * genome_one_filename;
   char               * genome_two_filename;
   char               * output_filename;
   char               * ref_exons_filename;
   char               * query_exons_filename;
   char               * ref_genes_filename;
   char               * query_genes_filename;
   char               * ref_gene_name;
   char               * query_gene_name;
   char               * ref_chrom;
   char               * query_chrom;
   double 		t, s, M;
   int 			T, x, p, u;
   unsigned int         l, v, Q, a, b, c, d;
   
 };

struct QGramOcc
 {
   unsigned int  occRef;
   unsigned int	 occQuery;
   unsigned int  length;
 };

struct MimOcc
 {
   unsigned int  startRef;
   unsigned int	 endRef;
   unsigned int  startQuery;
   unsigned int	 endQuery;
   unsigned int  error;
 };

struct PrevPos_L
 {
  unsigned int prev_L_ref;
  unsigned int prev_L_query;
 };

struct PrevPos_R
 {
  unsigned int prev_R_ref;
  unsigned int prev_R_query;
 };

typedef int32_t INT;

bool prefix(string str, string pref);
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
int remove_overlaps( vector<MimOcc> * mims, TSwitch sw );
double gettime ( void );
void usage ( void );
int alt_extend( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw, int alt );
int find_maximal_exact_matches( unsigned int l, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, TSwitch sw );
int editDistanceMyers( unsigned char * xInput, unsigned char * yInput );
int merge( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mims );
unsigned int rev_complement( unsigned char * str, unsigned char * str2, int iLen );
int adjust( unsigned int * edit_distance, int * q_start,  int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw );
int find_maximal_inexact_matches( TSwitch sw, unsigned char * ref, unsigned char * query, vector<QGramOcc> * q_grams, vector<MimOcc> * mnms, unsigned int qgram_size );
int extend( unsigned int * edit_distance,  int * q_start, int * q_end, int * r_start, int * r_end, unsigned char * xInput, unsigned char * yInput, TSwitch sw );
bool order(MimOcc a, MimOcc b);
unsigned int search( unsigned char * text, unsigned char * patt, unsigned int * score );
double scoring( MimOcc , unsigned char * ref, unsigned char * query );
