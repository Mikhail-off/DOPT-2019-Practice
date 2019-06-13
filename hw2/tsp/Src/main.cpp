#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include <set>
#include <fstream>
#include <iomanip> 
#include <random>

const size_t MAX_ITER = 70000000;
const size_t MEMORY_LIMIT = 7000;
const size_t EDGES_RATIO = 100;

typedef unsigned short int vertex_t;

struct Point {
	double x;
	double y;

	static double Distance( const Point& p1, const Point& p2 );
};

struct Edge {
	float Weight;
	vertex_t From;
	vertex_t To;
};

class CGeneticAlgo {
public:
	CGeneticAlgo( const std::vector<Point>& _points,
		const std::vector<std::vector<vertex_t>>& initialPopulation,
		double _selectionRatio, double _mutationRatio );
	void RunEvolution( std::vector<vertex_t>& hamPath );

	static double CalculatePathWeight( const std::vector<Point>& points, const std::vector<vertex_t>& path );
	
private:
	double selectionRatio;
	double mutationRatio;

	const std::vector<Point>& points;
	std::vector<std::vector<vertex_t>> population;

	std::mt19937 gen;
	std::uniform_real_distribution<double> randomGenerator;

	void selectBestSolution();
	void combineTwoPaths( const std::vector<vertex_t>& path1, const std::vector<vertex_t>& path2,
		std::vector<vertex_t>& hybrid );
	void mutatePath( std::vector<vertex_t>& path );
	bool pathComporator( const std::vector<vertex_t>& path1, const std::vector<vertex_t> path2 );
	bool canBeSolution( const std::vector<vertex_t>& path );

};

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Point
//////////////////////////////////////////////////////////////////////////////////////////////////

double Point::Distance( const Point& p1, const Point& p2 )
{
	return std::sqrt( ( p1.x - p2.x ) * ( p1.x - p2.x ) + ( p1.y - p2.y ) * ( p1.y - p2.y ) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Edge
//////////////////////////////////////////////////////////////////////////////////////////////////

bool operator < ( const Edge& e1, const Edge& e2 )
{
	return e1.Weight < e2.Weight;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void MST( const std::vector<Edge>& edges, std::vector<Edge>& treeEdges, size_t nVertex )
{
	std::vector<size_t> vertexId( nVertex );
	for( size_t i = 0; i < nVertex; ++i ) {
		vertexId[i] = i;
	}
	for( size_t i = 0; i < edges.size(); ++i ) {
		const Edge curEdge = edges[i];
		if( vertexId[curEdge.From] != vertexId[curEdge.To] ) {
			treeEdges.push_back( curEdge );
			size_t newId = vertexId[curEdge.From];
			size_t oldId = vertexId[curEdge.To];
			for( size_t j = 0; j < vertexId.size(); ++j ) {
				if( vertexId[j] == oldId ) {
					vertexId[j] = newId;
				}
			}
		}
		if( treeEdges.size() == nVertex - 1 ) {
			return;
		}
	}
}

std::vector<vertex_t> MakeHamiltonianPath( const std::vector<Edge>& mstEdges, const std::vector<Point>& points,
	size_t startVertex = 0 )
{
	std::vector<vertex_t> hamCycle;
	hamCycle.reserve( points.size() + 1 );
	size_t n = mstEdges.size() + 1;
	std::vector<bool> visited( n, false );

	hamCycle.push_back( startVertex );
	for( size_t i = 0; i < n; ++i ) {
		vertex_t last_vertex = hamCycle.back();
		visited[last_vertex] = true;
		bool foundNext = false;
		for( auto edge : mstEdges ) {
			if( edge.To == last_vertex ) std::swap( edge.To, edge.From );

			if( edge.From == last_vertex && visited[edge.To] == false ) {
				hamCycle.push_back( edge.To );
				foundNext = true;
				break;
			}
		}
		if( !foundNext ) {
			vertex_t nextVert = 0;
			double nextWeight = -1;
			for( vertex_t k = 0; k < points.size(); ++k ) {
				if( visited[k] ) continue;
				auto point = points[k];
				double curWeight = Point::Distance( point, points[last_vertex] );
				if( nextWeight == -1 || curWeight < nextWeight ) {
					nextWeight = curWeight;
					nextVert = k;
				}
			}
			hamCycle.push_back( nextVert );
		}
	}
	assert( hamCycle.size() == n + 1 );
	assert( hamCycle.back() == startVertex );
	hamCycle.pop_back();
	return hamCycle;
}

void MakeRandomReverse( const std::vector<Point>& points, std::vector<vertex_t>& hamPath )
{
	size_t n = hamPath.size();
	size_t left = std::rand() % n;
	if( left == 0 ) ++left;
	size_t right = left + std::rand() % ( n - left );
	if( left == right ) return;
	double diff = 0;
	diff += Point::Distance( points[hamPath[left - 1]], points[hamPath[left]] );
	diff += Point::Distance( points[hamPath[right]], points[hamPath[( right + 1 ) % n]] );
	diff -= Point::Distance( points[hamPath[left - 1]], points[hamPath[right]] );
	diff -= Point::Distance( points[hamPath[( right + 1 ) % n]], points[hamPath[left]] );
	if( diff > 1e-6 ) {
		std::reverse( hamPath.begin() + left, hamPath.begin() + right + 1 );
	}
}

void MakeRandomReverses( const std::vector<Point>& points, std::vector<vertex_t>& hamPath, int maxIter = MAX_ITER )
{
	assert( maxIter > 0 );
	for( int i = 0; i < maxIter; ++i ) {
		MakeRandomReverse( points, hamPath );
	}
}

void PrintAnswer( const std::vector<Point>& points, const std::vector<vertex_t>& path )
{
	std::cout << std::setprecision( 9 );
	std::cout << CGeneticAlgo::CalculatePathWeight( points, path ) << std::endl;
	for( auto vert : path ) {
		std::cout << vert + 1 << " ";
	}
	std::cout << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////// CGeneticAlgo
//////////////////////////////////////////////////////////////////////////////////////////////////

CGeneticAlgo::CGeneticAlgo( const std::vector<Point>& _points,
	const std::vector<std::vector<vertex_t>>& initialPopulation,
	double _selectionRatio, double _mutationRatio ) :
	points( _points ),
	population( initialPopulation.begin(), initialPopulation.end() ),
	selectionRatio( _selectionRatio ),
	mutationRatio( _mutationRatio ),
	randomGenerator( 0, 1 )
{
	std::random_device rd;
	gen = std::mt19937( rd() );
	/*
	for( size_t i = 0; i < initialPopulation.size(); ++i ) {
		size_t firstVertInd = 0;
		for( size_t j = 0; j < population[i].size(); ++i ) {
			if( population[i][j] == 0 ) {
				firstVertInd = j;
				break;
			}
		}
		std::rotate( population[i].begin(), population[i].end(), population[i].begin() + firstVertInd );
	}
	*/
}

double CGeneticAlgo::CalculatePathWeight( const std::vector<Point>& points, const std::vector<vertex_t>& path )
{
	if( path.size() <= 1 ) return 0;
	double weight = 0;
	for( size_t i = 1; i < path.size(); ++i ) {
		weight += Point::Distance( points[path[i]], points[path[i - 1]] );
	}
	weight += Point::Distance( points[path.front()], points[path.back()] );
	return weight;
}

bool CGeneticAlgo::canBeSolution( const std::vector<vertex_t>& path )
{
	return std::set<vertex_t>( path.begin(), path.end() ).size() == points.size();
}

bool CGeneticAlgo::pathComporator( const std::vector<vertex_t>& path1,
	const std::vector<vertex_t> path2 )
{
	double weight1 = CGeneticAlgo::CalculatePathWeight( points, path1 );
	double weight2 = CGeneticAlgo::CalculatePathWeight( points, path2 );
	return weight1 < weight2;
}

void CGeneticAlgo::RunEvolution( std::vector<vertex_t>& hamPath )
{
	std::sort( population.begin(), population.end(), pathComporator );
	size_t first = 0;
	size_t last = static_cast<size_t>( population.size() * selectionRatio );
	assert( first < last );
	assert( last < population.size() );
	for( size_t i = last; i < population.size(); ++i ) {
		size_t parent1 = std::rand() % last;
		size_t parent2 = parent1;
		while( parent1 == parent2 ) parent2 = std::rand() % last;
		combineTwoPaths( population[parent1], population[parent2], population[i] );
		mutatePath( population[i] );
	}
}

void CGeneticAlgo::combineTwoPaths( const std::vector<vertex_t>& path1, const std::vector<vertex_t>& path2,
	std::vector<vertex_t>& hybrid )
{
	size_t n = path1.size();
	assert( n == path2.size() );
	std::vector<bool> marked( n, false );
	size_t left = std::rand() % n;
	size_t right = left + std::rand() % ( n - left );
	size_t hybridIter = 0;
	for( size_t i = left; i < right; ++i ) {
		hybrid[hybridIter] = path1[i];
		marked[path1[i]] = true;
		hybridIter++;
	}
	for( size_t i = 0; i < n; ++i ) {
		if( !marked[path2[i]] ) {
			hybrid[hybridIter] = path2[i];
			hybridIter++;
		}
	}
	assert( hybridIter == n );
}


void CGeneticAlgo::mutatePath( std::vector<vertex_t>& path )
{
	if( randomGenerator( gen ) < mutationRatio ) {
		size_t left = std::rand() % path.size();
		size_t right = std::rand() % path.size();
		std::swap( path[left], path[right] );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
	// ��� �����, ����� ����� ���� � ���������� ���. ������ ���������� ���� �� ���� � �����
	std::ifstream input;
	std::ofstream output;
	if( argc > 1 ) {
		input = std::ifstream( argv[1] );
		std::cin.rdbuf( input.rdbuf() );
	}
	if( argc > 2 ) {
		output = std::ofstream( argv[2] );
		std::cout.rdbuf( output.rdbuf() );
	}


	size_t n = 0;
	std::cin >> n;
	std::vector<Edge> edges;
	std::vector<Point> points;
	points.reserve( n );
	edges.reserve( n * n / 2 / EDGES_RATIO );

	for( vertex_t i = 0; i < n; ++i ) {
		Point p;
		std::cin >> p.x >> p.y;

		for( vertex_t j = 0; j < points.size(); ++j ) {
			Edge edge;
			edge.From = j;
			edge.To = i;
			edge.Weight = static_cast<double>( Point::Distance( points[j], p ) );

			if( n > MEMORY_LIMIT ) {
				if( j % EDGES_RATIO == 0 ) {
					edges.push_back( edge );
				}
			} else {
				edges.push_back( edge );
			}

		}

		points.push_back( p );
	}
	//assert( false );
	std::sort( edges.begin(), edges.end() );
	std::vector<Edge> mstEdges;
	MST( edges, mstEdges, n );
	assert( mstEdges.size() == n - 1 );

	std::vector<vertex_t> hamPath = MakeHamiltonianPath( mstEdges, points, 0 );
	assert( std::set<vertex_t>( hamPath.begin(), hamPath.end() ).size() == n );

	// ���, ����!!! ����� ���������� ������� ������ �� �������
	std::vector<Edge>().swap( edges );

	//PrintAnswer(points, hamPath );
	MakeRandomReverses( points, hamPath );
	assert( std::set<size_t>( hamPath.begin(), hamPath.end() ).size() == n );

	//hamPath.push_back( hamPath.front() );
	PrintAnswer( points, hamPath );
	return 0;
}