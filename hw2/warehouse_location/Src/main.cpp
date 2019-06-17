#include "glpk.h"
#include <iostream>
#include <vector>
#include <fstream>


int main( int argc, char* argv[] )
{
	
	// Эта ересь, чтобы можно было в аргументах ком. строки передавать файл на вход и выход
	/*std::ifstream input;
	std::ofstream output;
	if( argc > 1 ) {
		input = std::ifstream( argv[1] );
		std::cin.rdbuf( input.rdbuf() );
	}
	if( argc > 2 ) {
		output = std::ofstream( argv[2] );
		std::cout.rdbuf( output.rdbuf() );
	}*/
	
	int n, m;
	std::cin >> n >> m;
	std::vector<double> capacity( n );
	std::vector<double> openCost( n );

	for( int i = 0; i < n; ++i ) {
		std::cin >> capacity[i] >> openCost[i];
	}

	std::vector<double> demand( m );
	for( int i = 0; i < m; ++i ) {
		std::cin >> demand[i];
	}
	std::vector<std::vector<double> > useCost( n, std::vector<double>( m ) );
	for( int i = 0; i < n; ++i ) {
		for( int j = 0; j < m; ++j ) {
			std::cin >> useCost[i][j];
		}
	}

	glp_prob* mip = glp_create_prob();
	glp_set_prob_name( mip, "   *_*   HARD PROBLEM!   *_*   " );
	glp_set_obj_dir( mip, GLP_MIN );

	// for all c : 
	//     sum_w x_wc = 1
	glp_add_rows( mip, n + m );
	for( int i = 1; i <= m; ++i ) {
		glp_set_row_bnds( mip, i, GLP_FX, 1.0, 1.0);
	}
	// for all w:
	//     sum_c x_wc * demand_c - y_w * capacyty_w <= 0
	for( int i = m + 1; i <= n + m; ++i ) {
		glp_set_row_bnds( mip, i, GLP_UP, 0.0, 0.0 );
	}

	// x_wc.size() = n * m
	// y_w.size() = n, y_w in {0, 1}
	glp_add_cols( mip, n * m + n );

	// for all w
	//     sum_c useCost * x_wc + openCost_w * y_w -> min
	for( int i = 1; i <= n; ++i ) {
		for( int j = 1; j <= m; ++j ) {
			int ind = ( i - 1 ) * m + j;
			glp_set_col_bnds( mip, ind, GLP_LO, 0.0, 0.0 );
			glp_set_obj_coef( mip, ind, useCost[i - 1][j - 1] );
		}
	}

	for( int i = 1; i <= n; ++i ) {
		int ind = m * n + i;
		glp_set_col_kind( mip, ind, GLP_BV );
		glp_set_obj_coef( mip, ind, openCost[i - 1] );
	}

	const int size = n * m + n * m + n;
	int ia[50000];
	int ja[50000];
	double ar[50000];

	// sum_w x_wc = 1
	int curVar = 1;
	for( int i = 1; i <= m; ++i ) {
		for( int j = 1; j <= n; ++j ) {
			ia[curVar] = i ;
			ja[curVar] = ( j - 1 ) * m + i;
			ar[curVar] = 1.0;
			curVar++;
		}
	}
	// for all w:
	//     sum_c x_wc * demand_c - y_w * capacyty_w <= 0
	for( int i = 1; i <= n; ++i ) {
		for( int j = 1; j <= m; ++j ) {
			ia[curVar] = m + i;
			ja[curVar] = ( i - 1 ) * m + j;
			ar[curVar] = demand[j - 1];
			curVar++;
		}
		// - y_w * capacity_w
		ia[curVar] = m + i;
		ja[curVar] = n * m + i;

		ar[curVar] = -capacity[i - 1];
		curVar++;
	}

	//std::cout << size << std::endl;
	glp_load_matrix( mip, size, ia, ja, ar );
	

	glp_iocp parm;
	glp_init_iocp( &parm );
	parm.presolve = GLP_ON;
	parm.msg_lev = GLP_MSG_OFF;
	parm.gmi_cuts = GLP_ON;
	parm.mir_cuts = GLP_ON;
	parm.mip_gap = 0.001;
	int err = glp_intopt( mip, &parm );

	std::vector<int> result;
	for( int i = 1; i <= n; ++i ) {
		if( glp_mip_col_val( mip, n * m + i ) ) {
			result.push_back( i );
		}
	}
	std::cout << result.size() << std::endl;
	for( int i = 0; i < result.size(); ++i ) {
		std::cout << result[i] << " ";
	}
	std::cout << std::endl;
	for( int i = 0; i < result.size(); ++i ) {
		for( int j = 1; j <= m; ++j ) {
			std::cout << glp_mip_col_val( mip, ( result[i] - 1 ) * m + j ) << " ";
		}
		std::cout << std::endl;
	}
	glp_delete_prob( mip );
	return 0;
}