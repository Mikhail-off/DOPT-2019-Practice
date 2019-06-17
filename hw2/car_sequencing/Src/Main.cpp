#include <iostream>
#include <vector>
#include <assert.h>


size_t ComputePenalty( const std::vector<size_t>& carLine,
	const std::vector<size_t>& robotWindow,
	const std::vector<size_t>& robotCapacity,
	const std::vector<std::vector<bool>> modelOptions )
{
	std::vector<size_t> penalties( robotCapacity.size(), 0 );
	std::vector<int> curReq( robotCapacity.size() );
	for( size_t option = 0; option < robotCapacity.size(); ++option ) {
		for( size_t i = 0; i < carLine.size(); ++i ) {
			if ( i + 1 )


			if( i + 1 >= robotWindow[option] ) {
				
			}
		}
	}
}

void GreedCarLine( const std::vector<size_t>& modelsCount,
	const std::vector<size_t>& robotCapacity,
	const std::vector<size_t>& robotWindow,
	const std::vector<std::vector<bool>>& modelOptions,
	std::vector<size_t>& carLine )
{
	size_t nCars = carLine.size();
	std::vector<size_t> carsRest( modelsCount.begin(), modelsCount.end() );
	std::vector<size_t> robotLoad( robotCapacity.size(), 0 );
	for( size_t i = 0; i < nCars; ++i ) {
		size_t modelToAdd = 0;
		bool canAddCar = true;
		size_t modelToAdd = 0;
		for( size_t model = 0; model < carsRest.size(); ++model ) {
			if( carsRest[model] == 0 ) continue;
			
			canAddCar = true;
			for( size_t option = 0; option < modelOptions[model].size(); ++option ) {
				if( robotCapacity[option] <= robotLoad[option] ) {
					canAddCar = false;
				}
			}

			if( canAddCar ) {
				modelToAdd = model;
				break;
			}
		}
		if( !canAddCar ) {
			while ( carsRest[modelToAdd++] == 0 ) {}
		}
		carLine[i] = modelToAdd;
		for( size_t option = 0; option < modelOptions[modelToAdd].size(); ++option ) {
			int toSub = 1;
			if( i + 1 < robotWindow[modelToAdd] ) {
				toSub = 0;
			} else {
				size_t modelToSub = carLine[i + 1 - robotWindow[modelToAdd]];
				toSub = modelOptions[modelToSub][option];
			}
			robotLoad[option] += modelOptions[modelToAdd][option];
			robotLoad[option] -= toSub;
		}
	}
}

int main()
{
	size_t nCars, nOptions, nModels;
	std::cin >> nCars >> nOptions >> nModels;

	std::vector<size_t> robotCapacity( nOptions );
	for( size_t i = 0; i < robotCapacity.size(); ++i ) {
		std::cin >> robotCapacity[i];
	}

	std::vector<size_t> robotWindow( nOptions );
	for( size_t i = 0; i < robotWindow.size(); ++i ) {
		std::cin >> robotWindow[i];
	}

	std::vector<size_t> modelsCount( nModels );
	std::vector<std::vector<bool>> modelsOptions( nModels,
		std::vector<bool>(nOptions));
	for( size_t i = 0; i < modelsCount.size(); ++i ) {
		std::cin >> modelsCount[i];
		for( size_t j = 0; j < nOptions; ++j ) {
			size_t option;
			std::cin >> option;
			modelsOptions[i][j] = static_cast<bool>( option );
		}
	}

	std::vector<size_t> carLine( nCars );
	GreedCarLine( modelsCount, robotCapacity, robotWindow, modelsOptions, carLine );
	
	return 0;
}