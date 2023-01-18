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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <iostream>
#include <float.h>
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>
#include "cnef.h"

static struct option long_options[] =
 {
   { "first-genome-file",              	required_argument, NULL, 'r' },
   { "sec-genome-file",              	required_argument, NULL, 'q' },
   { "output-file",             	required_argument, NULL, 'o' },
   { "min-seq-length",          	required_argument, NULL, 'l' },
   { "max-seq-length",          	required_argument, NULL, 'u' },
   { "sim-threshold",          		required_argument, NULL, 't' },
   { "ext-threshold",          		required_argument, NULL, 's' },
   { "threads", 			required_argument, NULL, 'T' },
   { "exons-file", 			required_argument, NULL, 'e' },
   { "ref-gene-file", 			optional_argument, NULL, 'g' },
   { "ref-gene-name", 			optional_argument, NULL, 'n' },
   { "ref-exons-file", 			required_argument, NULL, 'e' },
   { "query-gene-file", 		optional_argument, NULL, 'j' },
   { "query-gene-name", 		optional_argument, NULL, 'm' },
   { "query-exons-file", 		required_argument, NULL, 'f' },
   { "rev-complement",                  optional_argument, NULL, 'v' },
   { "remove-overlaps",			optional_argument, NULL, 'x' },
   { "ref-start", 			optional_argument, NULL, 'a' },
   { "ref-end", 			optional_argument, NULL, 'b' },
   { "query-start",                     optional_argument, NULL, 'c' },
   { "query-end",			optional_argument, NULL, 'd' },
   { "ref-chrom",                       optional_argument, NULL, 'y' },
   { "query-chrom",			optional_argument, NULL, 'z' },
   { "repeat-regions",			optional_argument, NULL, 'p' },
   { "merged-length",			optional_argument, NULL, 'M' },
   { "mem-length",			optional_argument, NULL, 'Q' },
   { "help",                    	no_argument,       NULL, 'h' },
   { NULL,                      	0,                 NULL,  0  }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> genome_one_filename            = NULL;
   sw -> genome_two_filename            = NULL;
   sw -> output_filename                = NULL;
   sw -> l                              = 50;
   sw -> u                              = 2000;
   sw -> t				= 1.0;
   sw -> s				= 0.05;
   sw -> v				= 0;
   sw -> x				= 1;
   sw -> ref_genes_filename		= NULL;
   sw -> ref_gene_name			= NULL;
   sw -> query_genes_filename		= NULL;
   sw -> query_gene_name		= NULL;
   sw -> ref_exons_filename		= NULL;
   sw -> query_exons_filename		= NULL;
   sw -> ref_chrom			= NULL;
   sw -> query_chrom			= NULL;
   sw -> a				= 0;
   sw -> b				= 0;
   sw -> c				= 0;
   sw -> d				= 0;
   sw -> p				= 1;
   sw -> T                              = 1;
   sw -> M				= 0.5;
   sw -> Q				= 18;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "q:r:o:e:f:g:j:x:n:m:l:u:t:s:v:a:b:c:d:y:z:p:T:M:Q:h", long_options, &oi ) ) != -1 ) 
    {

      switch ( opt )
       {

         case 'q':
           sw -> genome_two_filename= ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> genome_two_filename, optarg );
           args ++;
           break;

	 case 'r':
           sw -> genome_one_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> genome_one_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
          break;

	 case 'e':
           sw -> ref_exons_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> ref_exons_filename, optarg );
           args ++;
           break;

	case 'f':
           sw -> query_exons_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> query_exons_filename, optarg );
           args ++;
           break;

         case 'g':
           sw -> ref_genes_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> ref_genes_filename, optarg );
           args ++;
          break;

	  case 'j':
           sw -> query_genes_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> query_genes_filename, optarg );
           args ++;
          break;

         case 'n':
           sw -> ref_gene_name = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> ref_gene_name, optarg );
           args ++;
          break;

	  case 'm':
           sw ->  query_gene_name = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> query_gene_name, optarg );
           args ++;
          break;

         case 'y':
           sw -> ref_chrom = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> ref_chrom, optarg );
           args ++;
          break;

	  case 'z':
           sw ->  query_chrom = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> query_chrom, optarg );
           args ++;
          break;

         case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
           break;

	
         case 'u':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> u = val;
           break;

	  case 'x':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> x = val;
           break;

	case 't':
           val = atof( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> t = val;
           break;

	case 's':
           val = atof( optarg );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> s = val;
           break;

	case 'v':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> v = val;
           break;

	
	case 'a':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> a = val;
           break;

	
	case 'b':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> b = val;
           break;

	
	case 'c':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> c = val;
           break;

	
	case 'd':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> d = val;
           break;

	case 'p':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> p = val;
           break;
	
	 case 'T':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> T = val;
           break;

	 case 'M':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> M = val;
           break;

	 case 'Q':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> Q = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 7 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }

/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: CNEFinder <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -r, --ref-genome-file		<str>		FASTA reference genome filename.\n" );
   fprintf ( stdout, "  -q, --query-genome-file	<str>		FASTA query genome filename.\n" );
   fprintf ( stdout, "  -e, --exons-ref-file		<str>		GTF/GFF exon coordinates for reference genome filename.\n" );
   fprintf ( stdout, "  -f, --exons-query-file	<str>		GTF/GFF exon coordinates for query genome filename.\n" );
   fprintf ( stdout, "  -l, --min-seq-length		<int>		Minimum length of CNE.\n" );   
   fprintf ( stdout, "  -t, --sim-threshold		<dbl>		Threshold of similarity between sequences (0-1].\n" );
   fprintf ( stdout, "  -o, --output-file		<str>		Output filename with CNEs identified.\n\n" ); 

   fprintf ( stdout, "  Either 1. or 2.\n" );
   fprintf ( stdout, "    1.Search using gene name:\n" );
   
   fprintf ( stdout, "    -g, --ref-gene-file		<str>		GTF/GFF filename containing gene data for reference genome.\n" );
   fprintf ( stdout, "    -n, --ref-gene-name		<str>		Name of gene in reference genome in which CNEs will be identified.\n" );
   fprintf ( stdout, "    -j, --query-gene-file	<str>		GTF/GFF filename containing gene data for query genome.\n" );
   fprintf ( stdout, "    -m, --query-gene-name	<str>		Name of gene in query genome in which CNEs will be identified.\n\n" );
 
   fprintf ( stdout, "    2.Search using index position:\n" );
   fprintf ( stdout, "    -y, --ref-chrom		<str>		Chromosome of reference genome.\n" );
   fprintf ( stdout, "    -z, --query-chrom		<str>		Chromosome of query genome.\n" );
   fprintf ( stdout, "    -a, --ref-start		<int>		Start CNE search from this position of reference sequence.\n" );
   fprintf ( stdout, "    -b, --ref-end		<int>		End CNE search at this position of reference sequence.\n" );
   fprintf ( stdout, "    -c, --query-start		<int>		Start CNE search from this position of query sequence.\n" );
   fprintf ( stdout, "    -d, --query-end		<int>		End CNE search at this position of query sequence.\n\n" );

   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, "  -Q, --mem-length		<int>		Minimum length of maximal exact matches. Default:18.\n" );
   fprintf ( stdout, "  -M, --merged-length		<dbl>		Minimum length (in terms of CNE length) of merged matches to be extended. Default:0.5.\n" );
   fprintf ( stdout, "  -s, --ext-threshold		<dbl>		Threshold to further extend similarity threshold by. Default:0.05.\n" );
   fprintf ( stdout, "  -u, --max-seq-length		<int>		Set a maximum length for the CNE. Default:2000.\n" ); 
   fprintf ( stdout, "  -p, --repeat-regions		<int>		Choose 1 to filter repetitive regions of genomes or 0 otherwise. Default:1.\n");	
   fprintf ( stdout, "  -v, --rev-complement		<int>		Choose 1 to compute CNEs for reverse complement or 0 otherwise. Default:0.\n");						
   fprintf ( stdout, "  -x, --remove-overlaps		<int>		Choose 1 to remove overlapping CNEs or 0 otherwise. Default:1.\n\n" );  

   fprintf ( stdout, " Number of threads:\n" ); 
   fprintf ( stdout, "  -T, --threads			<int>		Number of threads to use. Default:1. \n" );
 }

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}
