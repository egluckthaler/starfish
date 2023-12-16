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
#include <vector>
#include <string.h>
#include <boost/algorithm/string.hpp>
#include <math.h>  
#include <sys/time.h>
#include <omp.h>
#include "cnef.h"

using namespace boost;

int main(int argc, char **argv)
{
	struct TSwitch  sw;

	FILE *          gen1_fd;                	 	
	FILE *          gen2_fd;    
	FILE *          ref_genes_fd;                	 	
	FILE *          ref_exons_fd;   
        FILE *          query_genes_fd;                	 	
	FILE *          query_exons_fd;       	 	
	FILE *          out_fd;                 	
        char *          genome_one_filename;          
	char *          genome_two_filename;         
	char *		ref_genes_filename;
	char * 		ref_exons_filename;
	char *		query_genes_filename;
	char * 		query_exons_filename;
        char *          output_filename;        
        char * 		ref_chrom;
	char * 		query_chrom;

        unsigned char ** genome1    = NULL;         	
	unsigned char ** genome2   = NULL;  

	unsigned char * ref   = NULL;  
	unsigned char * ref_id   = NULL;        	
	unsigned char * query   = NULL; 
	unsigned char * query_id   = NULL; 
     	
        unsigned char ** seq_id_genome1 = NULL;     	// the sequences id in memory
	unsigned char ** seq_id_genome2 = NULL;     	// the sequences id in memory

	unsigned char ** ref_genes = NULL;     	
	unsigned char ** ref_exons = NULL;   

	unsigned char ** query_genes = NULL;     	
	unsigned char ** query_exons = NULL;   

	unsigned int    i, j;
	unsigned int    q, l;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	omp_set_num_threads( sw.T );

	/* Check the arguments */
        if ( i < 7 )
        {
		usage ();
                return ( 1 );
        }
        else
        {

		genome_one_filename = sw . genome_one_filename;
		genome_two_filename = sw . genome_two_filename;

		if( sw . ref_genes_filename != NULL )
		{
			ref_genes_filename = sw . ref_genes_filename;	

			if ( ref_genes_filename == NULL  )
			{
				fprintf ( stderr, " Error: Cannot open gene data file for reference genome!\n" );
				return ( 1 );
			}
		}
		else if( sw . ref_chrom != NULL )
		{
			ref_chrom = sw . ref_chrom;	

			if ( ref_chrom == NULL )
			{
				fprintf ( stderr, " Error: Chromosome name for reference sequence required!\n" );
				return ( 1 );
			}
		}
		else
		{
			fprintf ( stderr, " Error: Choose gene name or index position for reference sequence to search for CNEs!\n" );
			return ( 1 );
		}

		
		if( sw . query_genes_filename != NULL )
		{
			query_genes_filename = sw . query_genes_filename;

			if ( query_genes_filename == NULL )
			{
				fprintf ( stderr, " Error: Cannot open gene data file for query genome!\n" );
				return ( 1 );
			}
		}
		else if( sw . query_chrom != NULL )
		{
			query_chrom = sw . query_chrom;	

			if ( query_chrom == NULL )
			{
				fprintf ( stderr, " Error: Chromosome name for query sequence required!\n" );
				return ( 1 );
			}
		}
		else
		{
			fprintf ( stderr, " Error: Choose gene name or index position for query sequence to search for CNEs!\n" );
			return ( 1 );
		}

		ref_exons_filename = sw . ref_exons_filename;
		query_exons_filename = sw . query_exons_filename;

		if ( genome_one_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for genome one!\n" );
			return ( 1 );
		}
		if ( genome_two_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for genome two!\n" );
			return ( 1 );
		}
		if ( ref_exons_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open exons data file for reference genome!\n" );
			return ( 1 );
		}
	
		if ( query_exons_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open exons data file for query genome!\n" );
			return ( 1 );
		}

		output_filename = sw . output_filename;
        }

	if( sw . l > sw . u )
	{
		fprintf ( stderr, " Error: Minimum length of CNE cannot be greater than maximum length of CNE!\n" );
		return ( 1 );
	}
	
	if( sw . a < 0 || sw . b < 0 || sw . c < 0 || sw . d < 0 )
	{
		fprintf ( stderr, "Error: Reference and Query positions cannot be less than 0" );
		return ( 1 );
	}

	/* Read the FASTA file for genome one in memory */
	fprintf ( stderr, " Reading the file: %s\n", genome_one_filename );
	if ( ! ( gen1_fd = fopen (  genome_one_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n",  genome_one_filename );
		return ( 1 );
	}

	char c;
        unsigned int num_seqs = 0;       
	unsigned int max_alloc_seq_id = 0;
	unsigned int max_alloc_seq = 0;
	c = fgetc( gen1_fd );

	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", genome_one_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs >= max_alloc_seq_id )
			{
				seq_id_genome1 = ( unsigned char ** ) realloc ( seq_id_genome1,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len = 0;
			unsigned int seq_id_len = 0;

			seq_id_genome1[ num_seqs ] = NULL;

			while ( ( c = fgetc( gen1_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id_len )
				{
					seq_id_genome1[ num_seqs ] = ( unsigned char * ) realloc ( seq_id_genome1[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len += ALLOC_SIZE;
				}
				seq_id_genome1[ num_seqs ][ seq_id_len++ ] = toupper( c );
			}
			seq_id_genome1[ num_seqs ][ seq_id_len ] = '\0';
			
		}

		if ( num_seqs >= max_alloc_seq )
		{
			genome1 = ( unsigned char ** ) realloc ( genome1,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq += ALLOC_SIZE;
		}

		unsigned int seq_len = 0;
		unsigned int max_alloc_seq_len = 0;

		genome1[ num_seqs ] = NULL;

		while ( ( c = fgetc( gen1_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n",  genome_one_filename );
				c = fgetc( gen1_fd );
				break;
			}
			if( c == '\n' || c == ' ' ) continue;

			if ( c == NA ) 
				c = '$';
			else if ( islower(c) && sw . p == 1 )
				c = '$';
			else c = toupper(c);

			if ( seq_len >= max_alloc_seq_len )
			{
				genome1[ num_seqs ] = ( unsigned char * ) realloc ( genome1[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}

			genome1[ num_seqs ][ seq_len++ ] = c;
		}

		if( seq_len != 0 )
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				genome1[ num_seqs ] = ( unsigned char * ) realloc ( genome1[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len += ALLOC_SIZE;
			}
			genome1[ num_seqs ][ seq_len ] = '\0';
			num_seqs++;	

		}
		
	} while( c != EOF );

	genome1[ num_seqs ] = NULL;

	if ( fclose ( gen1_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	/* Complete reading genome one */

		
	/* Read the FASTA file for genome two in memory */
	fprintf ( stderr, " Reading the file: %s\n", genome_two_filename );
	if ( ! ( gen2_fd = fopen ( genome_two_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", genome_two_filename );
		return ( 1 );
	}

	char cq;
        unsigned int num_seqs_q = 0;       
	unsigned int max_alloc_seq_id_q = 0;
	unsigned int max_alloc_seq_q = 0;

	cq = fgetc( gen2_fd );

	do
	{
		if ( cq != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", genome_two_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs_q >= max_alloc_seq_id_q )
			{
				seq_id_genome2 = ( unsigned char ** ) realloc ( seq_id_genome2,   ( max_alloc_seq_id_q + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id_q += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len_q = 0;
			unsigned int seq_id_len_q = 0;

			seq_id_genome2[ num_seqs_q ] = NULL;

			while ( ( cq = fgetc( gen2_fd ) ) != EOF && cq != '\n' )
			{
				if ( seq_id_len_q >= max_alloc_seq_id_len_q )
				{
					seq_id_genome2[ num_seqs_q ] = ( unsigned char * ) realloc ( seq_id_genome2[ num_seqs_q ],   ( max_alloc_seq_id_len_q + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len_q += ALLOC_SIZE;
				}
				seq_id_genome2[ num_seqs_q ][ seq_id_len_q++ ] =toupper( cq );
			}
			seq_id_genome2[ num_seqs_q ][ seq_id_len_q ] = '\0';
			
		}

		if ( num_seqs_q >= max_alloc_seq_q )
		{
			genome2 = ( unsigned char ** ) realloc ( genome2,   ( max_alloc_seq_q + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq_q += ALLOC_SIZE;
		}

		unsigned int max_alloc_seq_len_q = 0;
		unsigned int seq_len_q = 0;

		genome2[ num_seqs_q ] = NULL;

		while ( ( cq = fgetc( gen2_fd ) ) != EOF && cq != '>' )
		{
			if( seq_len_q == 0 && cq == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", genome_two_filename );
				cq = fgetc( gen2_fd );
				break;
			}
			if( cq == '\n' || cq == ' ' ) continue;

			if ( cq == NA ) 
				cq = '$';
			else if ( islower(cq) && sw . p == 1 )
				cq = '$';
			else cq = toupper(cq);

			if ( seq_len_q >= max_alloc_seq_len_q )
			{
				genome2[ num_seqs_q ] = ( unsigned char * ) realloc ( genome2[ num_seqs_q ],   ( max_alloc_seq_len_q + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len_q += ALLOC_SIZE;
			}

			genome2[ num_seqs_q ][ seq_len_q++ ] = cq;

		}

		if( seq_len_q != 0 )
		{
			if ( seq_len_q >= max_alloc_seq_len_q )
			{
				genome2[ num_seqs_q ] = ( unsigned char * ) realloc ( genome2[ num_seqs_q ],   ( max_alloc_seq_len_q + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len_q += ALLOC_SIZE;
			}
			genome2[ num_seqs_q ][ seq_len_q ] = '\0';
			num_seqs_q++;
		}
		
	} while( cq != EOF );

	genome2[ num_seqs_q ] = NULL;

	if ( fclose ( gen2_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	/* Complete reading genome two */


	/* Read the reference exons file */
	fprintf ( stderr, " Reading the file: %s\n", ref_exons_filename );
	if ( ! ( ref_exons_fd = fopen (  ref_exons_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n",  ref_exons_filename  );
		return ( 1 );
	}

	char ce;
        unsigned int num_seqs_e = 0;    
	unsigned int max_alloc_seq_e = 0;

	do
	{
		if ( num_seqs_e >= max_alloc_seq_e )
		{
			ref_exons = ( unsigned char ** ) realloc ( ref_exons,   ( max_alloc_seq_e + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq_e += ALLOC_SIZE;
		}

		unsigned int seq_len_e = 0;
		unsigned int max_alloc_seq_len_e = 0;

		ref_exons[ num_seqs_e ] = NULL;

		while ( ( ce = fgetc( ref_exons_fd ) ) != EOF )
		{
			if( seq_len_e == 0 && ce == '\n' )
			{
				fprintf ( stderr, " Omitting empty line in file %s!\n",  ref_exons_filename );
				ce = fgetc( ref_exons_fd );
				break;
			}

			if( ce == '\n' )
			{
				break;
			}
			else
			{

				if ( seq_len_e >= max_alloc_seq_len_e )
				{
					ref_exons[ num_seqs_e ] = ( unsigned char * ) realloc ( ref_exons[ num_seqs_e ],   ( max_alloc_seq_len_e + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_len_e += ALLOC_SIZE;
				}


				ref_exons[ num_seqs_e ][ seq_len_e++ ] = toupper( ce );
			}
	
		}

		if( seq_len_e != 0 )
		{
			if ( seq_len_e >= max_alloc_seq_len_e )
			{
				ref_exons[ num_seqs_e ] = ( unsigned char * ) realloc ( ref_exons[ num_seqs_e ],   ( max_alloc_seq_len_e + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len_e += ALLOC_SIZE;
			}
			ref_exons[ num_seqs_e ][ seq_len_e ] = '\0';
			num_seqs_e++;
		}
		
	} while( ce != EOF );


	ref_exons[ num_seqs_e ] = NULL;

	if ( fclose ( ref_exons_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	/* Complete reading referene exons file */


	unsigned int num_seqs_g = 0;   
	/* Read the reference genes file */
	if( sw . ref_gene_name != NULL )
	{
		fprintf ( stderr, " Reading the file: %s\n", ref_genes_filename );
		if ( ! ( ref_genes_fd = fopen (  ref_genes_filename, "r") ) )
		{
			fprintf ( stderr, " Error: Cannot open file %s!\n",  ref_genes_filename  );
			return ( 1 );
		}

		char cg;
		
		unsigned int max_alloc_seq_g = 0;

		do
		{
			if ( num_seqs_g >= max_alloc_seq_g )
			{
				ref_genes = ( unsigned char ** ) realloc ( ref_genes,   ( max_alloc_seq_g + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_g += ALLOC_SIZE;
			}

			unsigned int seq_len_g = 0;
			unsigned int max_alloc_seq_len_g = 0;

			ref_genes[ num_seqs_g ] = NULL;

			while ( ( cg = fgetc( ref_genes_fd ) ) != EOF )
			{
				if( seq_len_g == 0 && cg == '\n' )
				{
					fprintf ( stderr, " Omitting empty line in file %s!\n",  ref_genes_filename );
					cg = fgetc( ref_genes_fd );
					break;
				}

				if( cg == '\n' )
				{
					break;	
				}
				else
				{


					if ( seq_len_g >= max_alloc_seq_len_g )
					{
						ref_genes[ num_seqs_g ] = ( unsigned char * ) realloc ( ref_genes[ num_seqs_g ],   ( max_alloc_seq_len_g + ALLOC_SIZE ) * sizeof ( unsigned char ) );
						max_alloc_seq_len_g += ALLOC_SIZE;
					}

					ref_genes[ num_seqs_g ][ seq_len_g++ ] = toupper( cg );
				}
			}

			if( seq_len_g != 0 )
			{
				if ( seq_len_g >= max_alloc_seq_len_g )
				{
					ref_genes[ num_seqs_g ] = ( unsigned char * ) realloc ( ref_genes[ num_seqs_g ],   ( max_alloc_seq_len_g + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
					max_alloc_seq_len_g += ALLOC_SIZE;
				}
				ref_genes[ num_seqs_g ][ seq_len_g ] = '\0';

				num_seqs_g++;
			}
		
		} while( cg != EOF );


		ref_genes[ num_seqs_g ] = NULL;

		if ( fclose ( ref_genes_fd ) )
		{
			fprintf( stderr, " Error: file close error!\n");
			return ( 1 );
		}
	}
	/* Complete reading reference genes file */

	
	/* Read the query exons file */
	fprintf ( stderr, " Reading the file: %s\n", query_exons_filename );
	if ( ! ( query_exons_fd = fopen (  query_exons_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n",  query_exons_filename  );
		return ( 1 );
	}

	char cf;
        unsigned int num_seqs_f = 0; 
	unsigned int max_alloc_seq_f = 0;

	do
	{

		if ( num_seqs_f >= max_alloc_seq_f )
		{
			query_exons = ( unsigned char ** ) realloc ( query_exons,   ( max_alloc_seq_f + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq_f += ALLOC_SIZE;
		}

		unsigned int seq_len_f = 0;
		unsigned int max_alloc_seq_len_f = 0;

		query_exons[ num_seqs_f ] = NULL;

		while ( ( cf = fgetc( query_exons_fd ) ) != EOF )
		{
			if( seq_len_f == 0 && cf == '\n' )
			{
				fprintf ( stderr, " Omitting empty line in file %s!\n", query_exons_filename );
				cf = fgetc( query_exons_fd );
				break;
			}

			if( cf == '\n' )
			{
				break;
			}
			else
			{

				if ( seq_len_f >= max_alloc_seq_len_f )
				{
					query_exons[ num_seqs_f ] = ( unsigned char * ) realloc ( query_exons[ num_seqs_f ],   ( max_alloc_seq_len_f + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_len_f += ALLOC_SIZE;
				}


				query_exons[ num_seqs_f ][ seq_len_f++ ] = toupper ( cf );
			}
	
		}

		if( seq_len_f != 0 )
		{
			if ( seq_len_f >= max_alloc_seq_len_f )
			{
				query_exons[ num_seqs_f ] = ( unsigned char * ) realloc ( query_exons[ num_seqs_f ],   ( max_alloc_seq_len_f + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len_f += ALLOC_SIZE;
			}
			query_exons[ num_seqs_f ][ seq_len_f ] = '\0';

			num_seqs_f++;
		}
		
	} while( cf != EOF );


	query_exons[ num_seqs_f ] = NULL;

	if ( fclose ( query_exons_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	/* Complete reading query exons file */

	unsigned int num_seqs_j = 0; 
	/* Read the query genes file */
	if( sw . query_gene_name != NULL )
	{
		fprintf ( stderr, " Reading the file: %s\n", query_genes_filename );
		if ( ! ( query_genes_fd = fopen (  query_genes_filename, "r") ) )
		{
			fprintf ( stderr, " Error: Cannot open file %s!\n",  query_genes_filename  );
			return ( 1 );
		}

		char cj;
		      
		unsigned int max_alloc_seq_j = 0;
		do
		{	if ( num_seqs_j >= max_alloc_seq_j )
			{
				query_genes = ( unsigned char ** ) realloc (query_genes,   ( max_alloc_seq_j + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_j += ALLOC_SIZE;
			}

			unsigned int seq_len_j = 0;
			unsigned int max_alloc_seq_len_j = 0;

			query_genes[ num_seqs_j ] = NULL;

			while ( ( cj = fgetc( query_genes_fd ) ) != EOF )
			{
				if( seq_len_j == 0 && cj == '\n' )
				{
					fprintf ( stderr, " Omitting empty line in file %s!\n",  query_genes_filename );
					cj = fgetc( query_genes_fd );
					break;
				}

				if( cj == '\n' )
				{
					break;	
				}
				else
				{


					if ( seq_len_j >= max_alloc_seq_len_j )
					{
						query_genes[ num_seqs_j ] = ( unsigned char * ) realloc ( query_genes[ num_seqs_j ],   ( max_alloc_seq_len_j + ALLOC_SIZE ) * sizeof ( unsigned char ) );
						max_alloc_seq_len_j += ALLOC_SIZE;
					}

					query_genes[ num_seqs_j ][ seq_len_j++ ] = toupper( cj );
				}
			}

			if( seq_len_j != 0 )
			{
				if ( seq_len_j >= max_alloc_seq_len_j )
				{
					query_genes[ num_seqs_j ] = ( unsigned char * ) realloc ( query_genes[ num_seqs_j ],   ( max_alloc_seq_len_j + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
					max_alloc_seq_len_j += ALLOC_SIZE;
				}
				query_genes[ num_seqs_j ][ seq_len_j ] = '\0';
				num_seqs_j++;
			}
		
		} while( cj != EOF );


		query_genes[ num_seqs_j ] = NULL;

		if ( fclose ( query_genes_fd ) )
		{
			fprintf( stderr, " Error: file close error!\n");
			return ( 1 );
		}
	}
	/* Complete reading query genes file */


	fprintf ( stderr, " Pre-processing data\n" );


	if( sw . l < 4 )		
	{
		fprintf ( stderr, " Error: The minimum length of a match is 4.\n" );
		return ( 1 );	
	}


	if( sw . t < 0.7 )
	{
		fprintf( stderr, " Error: The similarity threshold cannot be less than 0.7!\n");
		return ( 1 );	
	}

	
	if( sw . t > 1 )
	{
		fprintf( stderr, " Error: The similarity threshold cannot be greater than 1.0!\n");
		return ( 1 );	
	}

	
	int geneOnePos = 0;
	int geneTwoPos = 0;

	unsigned int q_gram_size =  max( sw . Q/1.0, sw . l /  ( ( sw . l - floor( sw . t * sw . l ) ) + 1 ) ) ;
	sw . t = 1 - sw . t;

	if( sw . M <= 0 || sw . M > 1 )
	{
		fprintf( stderr, " Error: Minimum length of merged matches must be larger than 0 and smaller or equal to 1!\n");
		return ( 1 );
	}

	unsigned int start_genome_1 = 0;
	unsigned int end_genome_1 = 0;
	unsigned int start_genome_2 = 0;
	unsigned int end_genome_2 = 0;
	

	vector<int> * exons_g1_start = new vector<int>;
	vector<int> * exons_g1_end = new vector<int>;

	
	vector<int> * exons_g2_start = new vector<int>;
	vector<int> * exons_g2_end = new vector<int>;

	int countTabs = 0;
	string chromosome_g1;
	string chromosome_g1_nochr;
	string refGeneName = to_string( sw . a );
	refGeneName.append( " - ");
	refGeneName.append( to_string( sw . b ) );

	// Obtaining start and end coordinate from gene name

	if( sw . ref_gene_name != NULL )
	{
		refGeneName = reinterpret_cast<char*>( sw . ref_gene_name );
	
		for (i=0; i<refGeneName.length(); i++)
	    		refGeneName[i] = toupper( refGeneName[i] );

		refGeneName.append( "\t" );

		bool refGeneNameExists = false;
		for(int i=0; i<num_seqs_g; i++)
		{

			if( prefix( reinterpret_cast<char*>(ref_genes[i]), refGeneName ) == true )
			{
				geneOnePos  = i;
				refGeneNameExists = true;
				break;
			}	
		}

		if( refGeneNameExists == false )
		{
			fprintf( stderr, " Error: Reference gene name does not exist!\n");
			return ( 1 );
		}

		string startCoord_g1;
		
		string endCoord_g1;

		for(int i=0; i<strlen( (char * ) ref_genes[geneOnePos] ); i++)
		{
			string in = reinterpret_cast<char*>( ref_genes[geneOnePos] );

			if( ref_genes[geneOnePos][i] == '\t' )
				countTabs++;
			if( countTabs == 1 )
				chromosome_g1 +=  toupper(in[i]) ;
			if( countTabs == 2 )
				startCoord_g1 += in[i] ;
			if( countTabs == 3 )
				endCoord_g1 += in[i] ;
		}

		trim( startCoord_g1 );
		trim( endCoord_g1 );
		trim(  chromosome_g1 );
		start_genome_1 = atoi( startCoord_g1.c_str() );
		int diff = 0.05 * start_genome_1;
		start_genome_1 = start_genome_1 - diff;

		end_genome_1 = atoi( endCoord_g1.c_str() );
		end_genome_1 = end_genome_1 + diff;

		if( prefix( chromosome_g1 , "CHR" ) == false )
			chromosome_g1.insert( 0, "CHR" );
	
	}
	else
	{
		start_genome_1 = sw . a;
		end_genome_1 = sw . b;
		if( prefix( chromosome_g1 , "CHR" ) == false )
		{
			chromosome_g1_nochr.append( reinterpret_cast<char*>( sw . ref_chrom ) );
			chromosome_g1.insert( 0, "CHR" );
		}
		chromosome_g1.append( reinterpret_cast<char*>( sw . ref_chrom ) );
	}
	chromosome_g1_nochr.append( "\t" );
	chromosome_g1.append( "\t" );
	to_upper( chromosome_g1 );
	to_upper( chromosome_g1_nochr );

	string chromosome_g2;
	string chromosome_g2_nochr;
	string queryGeneName = to_string( sw . c );
	queryGeneName.append( " - ");
	queryGeneName.append( to_string( sw . d ) );
	if( sw . query_gene_name != NULL )
	{
		queryGeneName = reinterpret_cast<char*>(sw . query_gene_name );

		for (i=0; i<queryGeneName.length(); i++)
	    		queryGeneName[i] = toupper( queryGeneName[i] );
		queryGeneName.append( "\t" );

	
		bool queryGeneNameExists = false;
		for(int i=0; i<num_seqs_j; i++)
		{
		
			if( prefix(  reinterpret_cast<char*>(  query_genes[i] ) , queryGeneName ) == true )
			{
				geneTwoPos  = i;
				queryGeneNameExists = true;
				break;
			}
		}

		if( queryGeneNameExists == false )
		{
			fprintf( stderr, " Error: Query gene name does not exist!\n");
			return ( 1 );
		}


		countTabs = 0;
		string startCoord_g2;
		
		string endCoord_g2;
		for(int i=0; i<strlen( (char * ) query_genes[geneTwoPos] ); i++)
		{

			string in = reinterpret_cast<char*>(  query_genes[geneTwoPos] );
			if( query_genes[geneTwoPos][i] == '\t' )
				countTabs++;
			if( countTabs == 1 )
				chromosome_g2 += toupper(in[i]);
			if( countTabs == 2 )
				startCoord_g2 += in[i] ;
			if( countTabs == 3 )
				endCoord_g2 += in[i] ;
		}

		trim( startCoord_g2 );
		trim( endCoord_g2 );
		trim(  chromosome_g2 );
	

		start_genome_2 = atoi( startCoord_g2.c_str() );
		int diff = 0.05 * start_genome_2;
		start_genome_2 = start_genome_2 - diff;

		end_genome_2 = atoi( endCoord_g2.c_str() );
		end_genome_2 = end_genome_2 + diff;
		
		if( prefix( chromosome_g2, "CHR" ) == false )
			chromosome_g2.insert( 0, "CHR" );
	}
	else
	{
		start_genome_2 = sw . c;
		end_genome_2 = sw . d;
		if( prefix( chromosome_g2, "CHR" ) == false )
		{ 	chromosome_g2_nochr.append( reinterpret_cast<char*>( sw . query_chrom ) );
			chromosome_g2.insert( 0, "CHR" );

		}
		chromosome_g2.append( reinterpret_cast<char*>( sw . query_chrom ) );
	}
	chromosome_g2_nochr.append( "\t" );
	chromosome_g2.append( "\t" );
	to_upper( chromosome_g2 );
	to_upper( chromosome_g2_nochr );

	for ( i = 0; i < num_seqs_g; i ++ )
	{
		free ( ref_genes[i] );
	}
	for ( i = 0; i < num_seqs_j; i ++ )
	{
		free ( query_genes[i] );
	}
	free( ref_genes );
	free( query_genes );

	if( sw . a !=0 && start_genome_1 != sw . a )
	{
		fprintf( stderr, " Error: Start position for reference is different to that of gene\n. Search by either gene name or index position.\n  ");
		return ( 1 );
	}

	if( sw . b !=0 && end_genome_1 != sw . b )
	{
		fprintf( stderr, " Error: End position for reference is different to that of gene\n. Search by either gene name or index position.\n  ");
		return ( 1 );
	}

	
	if( sw . c !=0 && start_genome_2 != sw . c )
	{
		fprintf( stderr, " Error: Start position for query is different to that of gene\n. Search by either gene name or index position.\n"  );
		return ( 1 );
	}

	
	if( sw . d !=0 && end_genome_2 != sw . d )
	{
		fprintf( stderr, " Error: End position for query is different to that of gene\n. Search by either gene name or index position.\n " );
		return ( 1 );
	}


	if( sw . a >= sw . b && sw . a != 0 && sw . b !=0 )
	{
		fprintf( stderr, " Error: Start position of reference must be smaller than end position.\n " );
		return ( 1 );
	}
	if( sw . c >= sw . d && sw . c != 0 && sw . d !=0 )
	{
		fprintf( stderr, " Error: Start position of query must be smaller than end position.\n"  );
		return ( 1 );
	}
	
	//Obtaining exon coordinates
	

	for(int i=0; i<num_seqs_e; i++)
	{
		if( prefix( reinterpret_cast<char*>( ref_exons[i] ), chromosome_g1 ) == true ||  prefix( reinterpret_cast<char*>( ref_exons[i] ), chromosome_g1_nochr ) == true  )
		{
			countTabs = 0;
			string c1 = "";
			string c2 = "";

			string in = reinterpret_cast<char*>( ref_exons[i] );
			for(int j = 0; j<strlen( (char * ) ref_exons[i] ); j++)
			{
				if( ref_exons[i][j] == '\t' )
					countTabs++;
				if( countTabs == 1 )
					c1 += in[j];
				if( countTabs == 2 )
					c2 += in[j] ;
			}

			trim( c1 );
			trim( c2 );

			for(int k = 0; k<c1.length(); k++)
			{
				if( !isdigit( c1[k] ) ) 
				{	
					fprintf( stderr, " Error: Reference exon file format is incorrect.\n"  );
					return ( 1 );
				}
			}

			for(int k=0; k<c2.length(); k++ )	
			{
				if( !isdigit(  c2[k] ) ) 
				{
					fprintf( stderr, " Error: Reference exon file format is incorrect.\n"  );
					return ( 1 );
				}
			}

			if( start_genome_1 > 0 && end_genome_1 > 0 )
			{
				if( ( !( atoi( c2.c_str() ) < start_genome_1 )  && !(atoi( c1.c_str() ) > end_genome_1) ) )
				{
					exons_g1_start->push_back( atoi( c1.c_str() ) );
					exons_g1_end->push_back( atoi( c2.c_str() ) );
				}
			}
			else
			{
				exons_g1_start->push_back( atoi( c1.c_str() ) );
				exons_g1_end->push_back( atoi( c2.c_str() ) );

			}
		}

	}


	for(int i=0; i<num_seqs_f; i++)
	{
		if( prefix( reinterpret_cast<char*>(query_exons[i]), chromosome_g2 ) == true || prefix( reinterpret_cast<char*>( query_exons[i] ), chromosome_g2_nochr ) == true )
		{	
			countTabs = 0;
			string c1;
			string c2;
			string in = reinterpret_cast<char*>(  query_exons[i] );

			for(int j = 0; j<strlen( (char * ) query_exons[i] ); j++)
			{
				if( query_exons[i][j] == '\t' )
					countTabs++;
				if( countTabs == 1 )
					c1 += in[j] ;
				if( countTabs == 2 )
					c2 += in[j] ;
			}

			trim( c1 );
			trim( c2 );
	
			for(int k = 0; k<c1.length(); k++)
			{
				if( !isdigit( c1[k] ) )
				{
					fprintf( stderr, " Error: Query exon file format is incorrect.\n"  );
					return ( 1 );
				}
			}

			for(int k=0; k<c2.length(); k++ )	
			{
				if( !isdigit( c2[k] ) ) 
				{
					fprintf( stderr, " Error: Query exon file format is incorrect.\n"  );
					return ( 1 );
				}
			}

			if( start_genome_2 > 0 && end_genome_2 > 0 )
			{
				if( ( !( atoi( c2.c_str() ) < start_genome_2 )  && !(atoi( c1.c_str() ) > end_genome_2) ) )
				{	
					exons_g2_start->push_back( atoi( c1.c_str() ) );
					exons_g2_end->push_back( atoi( c2.c_str() ) );
				}
			}
			else
			{
				exons_g2_start->push_back( atoi( c1.c_str() ) );
				exons_g2_end->push_back( atoi( c2.c_str() ) );
			}
		}
	}

	//Obtain reference from genome1;
	
	bool g1Chromosome = false;
	for(int i=0; i<num_seqs; i++)
	{

		string chromosome = reinterpret_cast<char*>(seq_id_genome1[i]);
		chromosome.append( "\t" );
		
		if( chromosome == chromosome_g1 || chromosome == chromosome_g1_nochr )
		{

			if( sw . b == 0 )
				sw . b =  strlen( ( char* ) genome1[i] ) - 1;

			end_genome_1 = sw . b;

			ref = ( unsigned char * ) calloc ( ( end_genome_1 - start_genome_1 + 1 ) , sizeof( unsigned char ) );
			ref_id = ( unsigned char * ) calloc ( ( strlen( ( char* ) seq_id_genome1[i] ) + 1 ) , sizeof( unsigned char ) );

			if( sw . l > end_genome_1 - start_genome_1 )
			{
			
				fprintf( stderr, " Error: value of minimum length is larger than length of referene gene.\n" );
				return ( 1 );
			}

			if( start_genome_1 < 0 )
				start_genome_1 = 0;
		
			if( end_genome_1 > strlen( (char *) genome1[i] ) )
				end_genome_1 = strlen( (char *) genome1[i] ) - 1;
			
			if( start_genome_1 > strlen( (char*) genome1[i] ) || end_genome_1 > strlen( (char*) genome1[i] )  )
			{
			
				fprintf( stderr, " Error: Reference index position is larger than size of chromosome.\n" );
				return ( 1 );
			}

			memcpy( &ref[0], &genome1[i][start_genome_1], end_genome_1 - start_genome_1 );
			memcpy( &ref_id[0], &seq_id_genome1[i][0], strlen( ( char* ) seq_id_genome1[i] ));

			ref[ end_genome_1 - start_genome_1 ] = '\0';
			ref_id[ strlen( ( char* ) seq_id_genome1[i] ) ] = '\0';

			g1Chromosome = true;
			break;
		}
	}

	if( g1Chromosome == false )
	{
		fprintf( stderr, " Error: Chromosome %s not found in reference genome!\n", chromosome_g1.c_str() );
		return ( 1 );
	}


	//remove all exons from reference 
	string ref_str = reinterpret_cast<char*>(ref);
	#pragma omp parallel for
	for(int i=0; i<exons_g1_start->size(); i++)
		ref_str.replace(exons_g1_start->at(i) - start_genome_1, exons_g1_end->at(i) - exons_g1_start->at(i),  exons_g1_end->at(i) - exons_g1_start->at(i), '$' );
	
	ref =  (unsigned char*) ref_str.c_str();

	//Obtain query from genome2
	bool g2Chromosome = false;
	for(int i=0; i<num_seqs_q; i++)
	{

		string chromosome = reinterpret_cast<char*>(seq_id_genome2[i]);
		chromosome.append( "\t" );

		if( chromosome == chromosome_g2 || chromosome == chromosome_g2_nochr  )
		{
			if( sw . d == 0 )
				sw . d =  strlen( ( char* ) genome2[i] ) - 1;
			end_genome_2 = sw . d;
		
			query = ( unsigned char * ) calloc ( ( end_genome_2 - start_genome_2 + 1 ) , sizeof( unsigned char ) );
			query_id = ( unsigned char * ) calloc ( ( strlen( ( char* ) seq_id_genome2[i] ) + 1 ) , sizeof( unsigned char ) );

			if( sw . l > end_genome_2 - start_genome_2 )
			{
			
				fprintf( stderr, " Error: value of minimum length is larger than length of query gene.\n" );
				return ( 1 );
			}

			if( start_genome_2 < 0 )
				start_genome_2 = 0;
		
			if( end_genome_2 > strlen( (char *) genome2[i] ) )
				end_genome_2 = strlen( (char *) genome2[i] ) - 1;

			if( start_genome_2 > strlen( (char*) genome2[i] ) || end_genome_2 > strlen( (char*) genome2[i] )  )
			{
			
				fprintf( stderr, " Error: Query index position is larger than size of chromosome.\n" );
				return ( 1 );
			}

			memcpy( &query[0], &genome2[i][start_genome_2], end_genome_2 - start_genome_2);
			memcpy( &query_id[0], &seq_id_genome2[i][0], strlen( ( char* ) seq_id_genome2[i] ));

			query[ end_genome_2 - start_genome_2 ] = '\0';
			query_id[ strlen( ( char* ) seq_id_genome2[i] ) ] = '\0';

			g2Chromosome = true;
			break;
		}
	}

	if( g2Chromosome == false )
	{
		fprintf( stderr, " Error: Chromosome %s not found in reference genome!\n", chromosome_g2.c_str() );
		return ( 1 );
	}


	//removing exons from query
	string query_str = reinterpret_cast<char*>(query);

	#pragma omp parallel for
	for(int i=0; i<exons_g2_start->size(); i++)
		query_str.replace(exons_g2_start->at(i) - start_genome_2, exons_g2_end->at(i) - exons_g2_start->at(i),  exons_g2_end->at(i) - exons_g2_start->at(i), '$' );

	query = (unsigned char*) query_str.c_str();

	delete( exons_g1_start );
	delete( exons_g1_end );
	delete( exons_g2_start );
	delete( exons_g2_end );

	for ( i = 0; i < num_seqs; i ++ )
	{
		free ( genome1[i] );
		free( seq_id_genome1[i] );
	}
	for ( i = 0; i < num_seqs_q; i ++ )
	{
		free ( genome2[i] );
		free( seq_id_genome2[i] );
	}
	for ( i = 0; i < num_seqs_e; i ++ )
	{
		free ( ref_exons[i] );
	}
	for ( i = 0; i < num_seqs_f; i ++ )
	{
		free ( query_exons[i] );
	}
	free( genome1 );
	free( genome2 );
	free( seq_id_genome1 );
	free( seq_id_genome2 );
	free( ref_exons );
	free( query_exons );
	
	fprintf ( stderr, " Computing CNEs with minimum length %i, maximum length %i and similarity threshold %.2f\% \n", sw . l, sw . u, 100.0-sw.t * 100.0 );

	ofstream new_ref;
	new_ref.open(string(sw.output_filename)+"_new_ref.fa");
  	new_ref <<">"<<"new_ref_"+string(sw.output_filename)<<"\n"<<ref<<"\n";
  	new_ref.close();  

	ofstream new_query;
	new_query.open(string(sw.output_filename)+"_new_query.fa");
  	new_query <<">"<<"new_query_"+string(sw.output_filename)<<"\n"<<query<<"\n";
  	new_query.close();  

	vector<QGramOcc> * q_grams = new vector<QGramOcc>;
	vector<MimOcc> * mims = new vector<MimOcc>;

	double start = gettime();

	if( sw . v == 1 )
	{
		unsigned char * rc_seq = ( unsigned char * ) calloc ( ( strlen( ( char* ) query ) + 1 ) , sizeof( unsigned char ) );
			
		rev_complement( query, rc_seq , strlen( ( char* ) query ) );
		rc_seq[  strlen( ( char* ) query ) ] = '\0';

		find_maximal_exact_matches( q_gram_size , ref, rc_seq , q_grams, sw  );

		if( q_grams->size() == 0 && sw . l == 4  )
		{
			fprintf( stderr, " Error: No CNEs found.\n" );
			return ( 1 );
		}
		else if( q_grams->size() == 0  )
		{
			fprintf( stderr, " Error: No Matches found, try using a smaller value for minimum length.\n" );
			return ( 1 );
		}
		find_maximal_inexact_matches( sw , ref, rc_seq, q_grams, mims, q_gram_size );

		free( rc_seq );


	}
	else 
	{
		find_maximal_exact_matches( q_gram_size , ref, query , q_grams,  sw );

		if( q_grams->size() == 0 && sw . l == 4  )
		{
			fprintf( stderr, " Error: No CNEs found.\n" );
			return ( 1 );
		}
		else if( q_grams->size() == 0  )
		{
			fprintf( stderr, " Error: No Matches found, try using a smaller value for minimum length.\n" );
			return ( 1 );
		}

		find_maximal_inexact_matches( sw , ref, query, q_grams, mims, q_gram_size );
	}

	if( sw . x == 1 )
	{
		remove_overlaps( mims, sw );

	}

	delete( q_grams );

	fprintf ( stderr, " Preparing the output\n" );

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}

	
	for(int i=0; i<chromosome_g1.length(); i++)
		chromosome_g1[i] = tolower( chromosome_g1[i] );
	for(int i=0; i<chromosome_g2.length(); i++)
		chromosome_g2[i] = tolower( chromosome_g2[i] );

	trim( chromosome_g1 );
	trim( chromosome_g2 );

	//fprintf( out_fd, "%s%s%s%s%s%s%s\n", genome_one_filename, "\t", refGeneName.c_str(), "\t" , genome_two_filename,"\t", queryGeneName.c_str() );	
	for ( int i = 0; i < mims->size(); i++ )
	{
		if ( mims->at(i).endQuery - mims->at(i).startQuery >= sw . l || mims->at(i).endRef - mims->at(i).startRef >= sw . l )
		{

			unsigned int minLen = min( mims->at(i).endQuery - mims->at(i).startQuery, mims->at(i).endRef - mims->at(i).startRef);

			double threshold = (1.0 - (mims->at(i).error*1.0/minLen*1.0))*100.0;

			fprintf( out_fd, "%s%s%i%s%i%s%s%s%i%s%i%s%i%s%i%s%.2f\n", chromosome_g1.c_str(), "\t", mims->at(i).startRef+start_genome_1, "\t", mims->at(i).endRef + start_genome_1, "\t" , chromosome_g2.c_str() , "\t", mims->at(i).startQuery+start_genome_2, "\t", mims->at(i).endQuery+start_genome_2, "\t", mims->at(i).endRef - mims->at(i).startRef, "\t", mims->at(i).endQuery - mims->at(i).startQuery  ,"\t", threshold );
		}		
	}



	delete( mims );
		
	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}


	double end = gettime();

        fprintf( stderr, "Elapsed time: %lf secs.\n", end - start );
	ref_str.clear();
	free ( ref_id );
	query_str.clear();
	free ( query_id );
        free ( sw . genome_one_filename );
	free ( sw . genome_two_filename );
        free ( sw . output_filename );
	free ( sw . ref_exons_filename );
	free ( sw . query_exons_filename );
	free ( sw . ref_genes_filename );
	free ( sw . query_genes_filename );	
	free ( sw . ref_gene_name );
	free ( sw . query_gene_name );
	free ( sw . ref_chrom );
	free ( sw . query_chrom );

return 0;
}

unsigned int rev_complement( unsigned char * str, unsigned char * str2, int len )
{
	int i=0;
	while ( len >= 0 )
	{
		if ( str[len] == 'A' )
			str2[i++] = 'T';
		else if( str[len] == 'C')
			str2[i++] = 'G';
		else if( str[len] == 'G')
			str2[i++] = 'C';
		else if( str[len] == 'T')
			str2[i++] = 'A';
		else if( str[len] == 'N')
			str2[i++] = 'N';
		len--;
	}
return 1;
}


bool prefix(string str, string pref)
{
	if( strncmp(str.c_str(), pref.c_str(), pref.length() ) == 0 )
		return true;

return false;
}
