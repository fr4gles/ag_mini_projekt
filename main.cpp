#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ga/ga.h>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <ctime>

#define MAX(x,y) (x)>(y)?(x):(y)
#define MIN(x,y) (x)<(y)?(x):(y)

class Block;

int popsize								= 1000; // wielkosc populacji
int ngen								= 3000; // ilosc generacji
float pmut								= 0.05; // prawdopodobienstwo mutacji
float pcross							= 0.65; // prawdopodobienstwo krzyzowania

int kolejnosc							= 0;	// zmienna uzywana przy inicjalizacji osobnika - 1 osobnik jest inicjalizowany inaczej							
int przewrocil_sie_global				= 0;	// zmienna okreslajaca max w generacjach liczbe klockow po jaka zostala dobrze ulozona 
int przewrocil_sie_global_w_generacji	= 0;	// zmienna okreslajaca max w danej generacji liczbe klockow po jaka zostala dobrze ulozona 
//int max_width = 0;

unsigned int amount_of_blocks			= 0;	// liczba wczytanych klocków
std::vector<Block> blocks_from_file;			// tablica wczytancyh klocków

/// klasa Block -  klasa klocków
/// - s³u¿y do przechowywania wartosci klocków 
/// - przechowuje wszystkie potrzebne wartosci do okreslenia po³o¿enia, wielkoœci oraz obrotu klocka
class Block
{
public:
	/// domyslny konstruktor, inicjalizacja danych
	Block()
	{
		i_k		= 0;
		w_k		= 0;
		h_k		= 0;
		x_k		= 0.0;
		turned	= false;
		mid_point = 0.0;
		setHalfPoint();
	}

	/// konsktruktor z ustawieniem wartosci - wiekszosc z nich moze byc domyslna
	Block(const int &i, const int &w = 0, const int &h = 0, const float &x = 0, const bool& t = 0, const float &mid=0.0)
	{
		i_k		= i;
		w_k		= w;
		h_k		= h;
		x_k		= x;
		turned	= t; 
		mid_point = mid;
		setHalfPoint();
	}

	/// fukncja ustawiajaca wartosci klockow - odpowiednik konstruktora
	void setBlock(const int &i, const int &w, const int &h, const float &x, const bool& t, const float &mid=0.0)
	{
		i_k		= i;
		w_k		= w;
		h_k		= h;
		x_k		= x;
		turned	= t;
		mid_point = mid;
		setHalfPoint();
	}

	/// funkcja kopiujaca wartosci z zdanego klocka
	/// uzywana przy inicjalizacji osobnika
	void cloneBlock(const Block &block)
	{
		i_k		= block.i_k;
		w_k		= block.w_k;
		h_k		= block.h_k;
		x_k		= block.x_k;
		mid_point = block.mid_point;
		turned	= block.turned;
		setHalfPoint();
	}

	/// funkcja obliczajaca polowe podstawy
	/// potrzeba do obliczen srodka ciezkosci klocka
	void setHalfPoint() const
	{
		half_w	= w_k/2.0;
	}

	/// funckja obracajaca klocek
	void rotate() const 
	{
		int tmp = 0;
		tmp = this->h_k;
		this->h_k = this->w_k;
		this->w_k = tmp;

		if(turned)
			this->turned = false;
		else
			this->turned = true;
		setHalfPoint();
	}

//////////////////////////////////////////////////////////////////////////
	///operatory na potrzeby biblioteki
	operator int()
	{
		return this->i_k;
	}

	int operator!=(const Block& tmp) const
	{
		return this->i_k != tmp.i_k;
	}

	int operator==(const Block& tmp) const
	{
		return this->i_k == tmp.i_k;
	}

	friend std::ostream& operator<< (std::ostream &out, const Block& block);
//////////////////////////////////////////////////////////////////////////

	/// deklaracja zmiennych 
	mutable int		i_k,		/// indeks klocka
					w_k,		/// szerokosc klocka
					h_k;		/// wysokosc klocka
	mutable float	x_k,		/// przesuniecie klocka wzgledem poprzedniej lewej krawedzi
					mid_point,	/// wspolrzedna srodka klocka ( na osi X )
					half_w;		/// polowa szerokosci klocka
	mutable bool	turned;		/// czy klocek zostal obrocony 
};

std::ostream& operator<<(std::ostream &out, const Block& block)
{
	return (out << block.i_k);
}


//////////////////////////////////////////////////////////////////////////
/////////////////////// deklaracje funkcji ///////////////////////////////
//////////////////////////////////////////////////////////////////////////
/// funkcja dostosowania - dalszy opis w ciele funkcji
float objective(GAGenome &);

/// funkcja inicjalizujaca osobnika - osobnik jest inicjalizowany klockami w losowej kolejnosci
void init_my_population(GAGenome &ga);

/// funkcja na potrzeby testow sluzaca do generowana pliktu testowego z klockami
/// - nazwa pliku, ilosc klockow
void generate_blocks_file(const std::string &, const int &);

/// funckja zczytujaca klocki z pliku 
/// - nazwa pliku
void read_blocks_file(const std::string &);

/// funckja generujaca losowa liczbe zmiennoprzecinkowa typu float
/// - przedzial od, do
float random_float(const float &, const float &);

/// funkcja zapisujaca wyniki do pliku
void write_blocks_file(const std::string &, const GA1DArrayGenome<Block>&);

/// funkcja napisana na potrzeby inicjalizacji - jeden z osobnikow jest inicjalizawany posortowanymi wartosciami
/// wprowadza roznorodnosc w populacji
bool sortuj_dobrze (const Block &i,const Block &j) { return (i.w_k>j.w_k); }

/// funkcja obliczajaca szerokosc najlepszego osobnika
/// - na potrzeby testow
float maksymalna_szerokosc_osobnika(const GA1DArrayGenome<Block>& );
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/////////////////////////// main /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	std::time_t start = clock();	/// pobranie czasu startu programu
	srand((unsigned)time(0));		/// wymuszenie losowania roznych wartosci w rand()	

	//generate_blocks_file("plik.txt", 300);/// generowanie pliku z klockami

	read_blocks_file("plik.txt");

//////////////////////////////////////////////////////////////////////////
	/// osobnik
	/// -  kazdy osobnik sklada sie z wczytanej ilosci klockow
	/// - jest inicjalizowany losowym ulozeniem klockow
	///		- krzyzowanie	- OrderCrossOver
	///		- mutator		- SwapMutator
	GA1DArrayGenome<Block> genome(amount_of_blocks,objective);	
	genome.initializer(init_my_population);
	genome.crossover(GA1DArrayGenome<Block>::OrderCrossover);
	genome.mutator(GA1DArrayGenome<Block>::SwapMutator);
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
	/// algorytm SGA - simple genetic algorithm
	/// metoda selekcji - RANKSELECTOR
	/// prawdopodobienstwa mutacji oraz krzyzowania zgodnie z ustawieniami
	/// zastosowano skalowanie wynikow funkcji objective
	GASimpleGA ga(genome);
	GARankSelector select;

	ga.selector(select);
	ga.populationSize(popsize);
	ga.nGenerations(ngen);
	ga.pMutation(pmut);
	ga.pCrossover(pcross);

	//GANoScaling scaling;
	GASigmaTruncationScaling scaling;
	ga.scaling(scaling);

	ga.scoreFrequency(1);
	ga.flushFrequency(100);
	ga.initialize((unsigned)time(0));
//////////////////////////////////////////////////////////////////////////
	
	int tmp_ilosc_powtorzen = 0;			/// zmienna tmp, okresla ilosc pwotorzen maksymalnej ilosci polozonych klockow klockow
	int tmp_przewrocil_sie_global = 0;		/// zmienna tmp, okresla ilosc polozonych klockow w poprzedniej generacji
	int wykryto_stagnacje = 0;				/// zmienna okreslajaca stopien stagnacji - nic sie nie zmienia, nie mozna ulozyc wiecej klockow
	int ktora_gen = 0;						/// zmienna okresla numer generacji

//////////////////////////////////////////////////////////////////////////
	/// petla while - glowny przebieg programu
	/// rozwiazywanie problemu - serce algorytmu genetycznego
	while(!ga.done())
	{
	
	//////////////////////////////////////////////////////////////////////////
		/// ta czesc sluzy do okreslenia jak wiele razy algorytm nie byl w stanie znalezc lepszego wyniku
		if(przewrocil_sie_global==tmp_przewrocil_sie_global)
			++tmp_ilosc_powtorzen;
		else
			tmp_ilosc_powtorzen = 0;

		tmp_przewrocil_sie_global = przewrocil_sie_global;
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
		ga.maximize();		/// maksymalizowanie wynikow
		ga.step();			/// nastepna generacja ...
	//////////////////////////////////////////////////////////////////////////

	
	//////////////////////////////////////////////////////////////////////////
		/// wypis obecnych wynikow
		std::cout << ktora_gen << " :   " << przewrocil_sie_global_w_generacji << "  :  " << przewrocil_sie_global << "      \r" ;
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
		/// czesc odpowiadajaca za zmiane pradopodobienstw krzyowania oraz mutacji
		/// w wypadku gdy algorytm wiele razy nie znalazl dobrego rozwiazania
		/// - w zalozeniu ma pomoc znalezc lepsze rozwiazanie oraz zapobiegac stagnacji w populacji
		///		- raz wartosci prawdopodobienstw sa zwiekszane, raz sa liczba losowa
		if(tmp_ilosc_powtorzen>20)
		{
			++wykryto_stagnacje;
			if(pmut < 0.09)
				pmut += 0.01;
			else
				pmut = random_float(0.5,1.0);

			if(pcross < 1.0)
				pcross += 0.05;
			else
				pcross = random_float(0.0,0.5);

			ga.pMutation(pmut);
			ga.pCrossover(pcross);
			
			tmp_ilosc_powtorzen = 0;
		}

		if(wykryto_stagnacje > 10)
		{
			pmut = 0.0;
			pcross = 0.0;
			ga.pMutation(pmut);
			ga.pCrossover(pcross);
			wykryto_stagnacje = 0;
		}
	//////////////////////////////////////////////////////////////////////////
		


	//////////////////////////////////////////////////////////////////////////
		/// ustawianie wartosci zmiennych 
		przewrocil_sie_global_w_generacji = 0;
		ktora_gen++;
	//////////////////////////////////////////////////////////////////////////
	}

//////////////////////////////////////////////////////////////////////////
	/// zapis najleszego osobnika do pliku
	GA1DArrayGenome<Block>& best_of_all = (GA1DArrayGenome<Block> &)ga.statistics().bestIndividual();
	write_blocks_file("out.txt", best_of_all/*(GA1DArrayGenome<Block> &)ga.statistics().bestPopulation().worst()*/);
	std::cout << "\n" << maksymalna_szerokosc_osobnika(best_of_all) << "     " << best_of_all.score() <<" \n";
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
	/// konczenie programu wraz z wysietleniem czasu wykonia
	std::time_t stop = clock();
	std::cout << "Program wykonal sie w :" << ((float)(stop-start)/CLOCKS_PER_SEC) << "[s]" << "\n";
//////////////////////////////////////////////////////////////////////////

	return 0;
}

float objective(GAGenome & c)
{
	GA1DArrayGenome<Block> & genome = (GA1DArrayGenome<Block> &)c;

	float	result = 0.0,
			punkt_srodka_ciezkosci = 0.0,
			srodek_ciezkosci = 0.0,
			max_prawe_wychylenie = 0.0,
			x1 = 0.0,
			x2 = 0.0,
			y_k,
			kara = 0.0,
			kara_x = random_float(0.0,0.5),
			kara_y = random_float(0.0,0.5),
			rzut_moneta = 0.0;

	int		przewrocil_sie = 0,
			rotated_blocks_good = 0,
			rotated_blocks_bad	= 0,
			local_max_width = 0;
	
	for(unsigned int i=0; i<amount_of_blocks;++i)
	{
		///obrót klocków - 50% szans na obrót
		rzut_moneta = random_float(0.0,1.0);

		if(rzut_moneta >= 0.5)
		{
			//if(genome.gene(i).w_k < genome.gene(i).h_k)
				genome.gene(i).rotate();
		}

		if(genome.gene(i).w_k > genome.gene(i).h_k)
			++rotated_blocks_good;
		else
			++rotated_blocks_bad;
		//////////////////////////////////////////////////////////////////////////
		switch(i)
		{
			case 0: 
				genome.gene(i).x_k = 0.0;
				x1 = genome.gene(i).half_w;
				x2 = -x1;
				genome.gene(i).mid_point = 0.0;
			break;
			default: 
				/// wyznaczanie x_k
				max_prawe_wychylenie = genome.gene(i-1).w_k - genome.gene(i).half_w;
				genome.gene(i).x_k = random_float((-1.0)*genome.gene(i).half_w+(kara_x*genome.gene(i).half_w), max_prawe_wychylenie-(kara_y*max_prawe_wychylenie));

				/// wyznaczanie wspolrzedneych srodkowego punktu
				y_k = genome.gene(i).w_k - ( (-1.0)*(genome.gene(i).x_k) + genome.gene(i-1).w_k );
				genome.gene(i).mid_point = genome.gene(i-1).mid_point + ((genome.gene(i).x_k + y_k)/2.0);	
			break;
		}

		if(kolejnosc)
			local_max_width =  MAX(local_max_width,genome.gene(i).w_k);
	}

	int n=0;
	for(int i = 1;i<amount_of_blocks;++i)
	{
		if(przewrocil_sie)
			break;
		n=0; 
		srodek_ciezkosci = 0.0;
		punkt_srodka_ciezkosci = 0.0;

		switch(i)
		{
		case 0 :	if(genome.gene(i).w_k<genome.gene(i).h_k)
						kara -= amount_of_blocks;
					
					if(genome.gene(i).w_k != local_max_width)
						kara -= 2.0*amount_of_blocks;
			break;
		default:	if(genome.gene(i).w_k<genome.gene(i).h_k)
						kara -= (amount_of_blocks-i);
			break;
		}

		for (int j = i; j>0;--j)
		{
			if(przewrocil_sie)
				break;

			/// obliczanie szerokosci wiezy
			x1 = MIN(genome.gene(i).mid_point - genome.gene(i).half_w,x1);
			x2 = MAX(genome.gene(i).mid_point + genome.gene(i).half_w,x2); 

			/*std::cout << x1 << "    " << x2 << "   " << x2+fabs(x1) << std::endl;*/
			//////////////////////////////////////////////////////////////////////////
			++n;
			punkt_srodka_ciezkosci += genome.gene(j).mid_point;
			srodek_ciezkosci = (punkt_srodka_ciezkosci) / n;

			switch(przewrocil_sie)
			{
			case 0:
				if( ( genome.gene(j).x_k >= 0 ) && ( (srodek_ciezkosci > ((genome.gene(j-1).mid_point + genome.gene(j-1).half_w))) || (srodek_ciezkosci < ((genome.gene(j-1).mid_point - genome.gene(j-1).half_w))) ))
				{
					przewrocil_sie = j;
					genome.gene(j).x_k = (float)INT_MAX;
				}
				else if( ( genome.gene(j).x_k < 0 ) && ( (srodek_ciezkosci > ((genome.gene(j-1).mid_point + genome.gene(j-1).half_w))) || (srodek_ciezkosci < ((genome.gene(j-1).mid_point - genome.gene(j-1).half_w))) ))
				{
					przewrocil_sie = j;
					genome.gene(j).x_k = (float)INT_MAX;
				}
				break;
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	przewrocil_sie_global = MAX(przewrocil_sie,przewrocil_sie_global);
	przewrocil_sie_global_w_generacji = MAX(przewrocil_sie,przewrocil_sie_global_w_generacji);

	
	//	result = 10.0*amount_of_blocks + 0.5*(x2+fabs(x1)) + 2.0*(rotated_blocks_good-rotated_blocks_bad) +	0.01*kara;
	if(przewrocil_sie==0)
		przewrocil_sie = amount_of_blocks;

	float wspolczynnik = 0.5;

	if(amount_of_blocks<100)
		wspolczynnik = 8.0;

	result	=	70.5*przewrocil_sie
			+	wspolczynnik*(x2+fabs(x1));

	if(amount_of_blocks > 50)
	{
		result	+= 	rotated_blocks_good
				-	2.0*rotated_blocks_bad
				+	0.5*kara;
	}

	//std::cout << x1 << "    " << x2 << "   " << x2+fabs(x1) << "    " << result <<std::endl;

	return result;
}

float maksymalna_szerokosc_osobnika(const GA1DArrayGenome<Block>& best)
{
	float	x1 = 0.0,
			x2 = 0.0;

	for(unsigned i=0; i<amount_of_blocks; ++i)
	{
		if(best.gene(i).x_k ==(float)INT_MAX) 
			break;

		x1 = MIN(best.gene(i).mid_point - best.gene(i).half_w,x1);
		x2 = MAX(best.gene(i).mid_point + best.gene(i).half_w,x2); 
	}

	return x2+fabs(x1);
}

void init_my_population(GAGenome &ga)
{
	GA1DArrayGenome<Block> &my_gene = (GA1DArrayGenome<Block> &)ga;

	if(kolejnosc)
		std::random_shuffle(blocks_from_file.begin(),blocks_from_file.end());

	int tmp = 0;
	for(unsigned int i=0;i<amount_of_blocks;++i)
	{
		my_gene[i].cloneBlock(blocks_from_file[i]);
	}
	
	kolejnosc++;
}

void generate_blocks_file(const std::string &nazwa, const int &N)
{
	std::ofstream blocks_file(nazwa.c_str());

	blocks_file << N << std::endl;

	for(unsigned i=0; i<N; ++i)
	{
		blocks_file << i << " " << rand()%16+5 << " " << rand()%16+5 << std::endl;
	}
}

void read_blocks_file(const std::string &tmp) 
{
	blocks_from_file.reserve(301);
	FILE *file;
	file = fopen(tmp.c_str(),"r");
	if(file==NULL)	{ printf("Blad otwierania pliku z klockami. ABORT\n\n"); exit(-1); }

	int i,w,h; 
	fscanf(file,"%d",&amount_of_blocks);
	for(unsigned int j=0;j<amount_of_blocks;++j)
	{
		fscanf(file,"%d %d %d",&i,&w,&h);
		blocks_from_file.push_back(Block(i,w,h));
	}

	std::sort(blocks_from_file.begin(), blocks_from_file.end(),sortuj_dobrze);
}

float random_float(const float &LO, const float &HI)
{
	return LO + ((float)rand()/((float)RAND_MAX/(HI-LO)));
}

void write_blocks_file(const std::string &name, const GA1DArrayGenome<Block>& data)
{
	std::ofstream file(name.c_str());

	for(unsigned j=0;j<amount_of_blocks;++j)
	{
		file << data.gene(j).i_k << "\t" << data.gene(j).turned << "\t" << data.gene(j).x_k << std::endl;
	}

	file.close();
}

