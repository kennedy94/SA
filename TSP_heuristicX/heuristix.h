#include <iostream>
#include <fstream>
#include <limits>
#include <math.h>
using namespace std;



class TSP_instance{
private:
	double **D;	//Matriz do problema;
	int n;		//#cidades
	int *sopt;	//solucao otima ate o momento
	int salto;
	int *set; //vetor de binarios, determina se 'i' ja esta na solucao
	int *sol_aleatoria(int s[]);
public:
	void opt4(int ss[]);
	double f(int ss[]);
	void SimulatedA();
	double custo_opt;
	// OBS: maD Define se D é triang. superior (1), inferior (0) ou Cheia (2).
	//construtor para ler a instancia
	TSP_instance(const char* filename, int matD){

		ifstream input(filename, ifstream::in);
        if (input.fail()) {
            cerr << "		Arquivo \"" << filename << "\" nao encontrado." << endl;
			getchar();
            exit(1);
        }else if (matD < 0 || matD > 3){
			cerr << "		Defina um tipo de matriz válido" << endl;
			getchar();
			exit(2);
		}

		int i, j;

		input >> n;

		salto = n/50; //definindo salto
		iniciar_variaveis(); //iniciando variaveis

		cout << "Problema com n = " << n << " e salto = " << salto << endl; 

		// Entrada de dados para matriz triangular superior
		switch(matD){
		case 0:
			for (i = 1; i < n; i++){
				for (j = 0; j <= i-1; j++){
					input >> D[i][j];
					D[j][i]= D[i][j];
				}
			}
			break;
		case 1:
			for (i = 0; i < n-1; i++){
				for (j = i+1; j < n; j++){
					input >> D[i][j];
					D[j][i]= D[i][j];
				}
			}
			break;
		case 2:
			for (i = 0; i < n; i++){
				for (j = 0; j < n; j++){
					input >> D[i][j];
				}
			}
			break;
		case 3://OBS: ESSA EH A MATRIZ GPS
			double *x = new double[n];
			double *y = new double[n];
			for (i = 0; i < n; i++) {
				input>> j;
				input >> x[i];
				input >> y[i];
				x[i]--;
				y[i]--;
			}
			for (i = 0; i < n; i++){
				for (j = 0; j < n; j++){
					if (i!=j) 
						D[i][j]= sqrt( pow(x[i] - x[j],2.0) + pow(y[i] - y[j], 2.0) );
				}
			}
			delete x,y;
			break;
		}

		
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				if (i==j) D[i][j]= 200000;
				//cout << D[i][j] << " ";
			}
			//cout << endl;
		}
		
		input.close();
	}

	void TSP_instance::HeuristicaRUIM();
	void iniciar_variaveis(){
		custo_opt = INT_MAX;
		D = new double*[n];
		for(int i = 0; i < n; i++)
			D[i] = new double[n];
		sopt = new int[n];
		set = new int[n];
	}

	void g(int ss[])
	{
		int i;
		double custo=0;
		//for (i = 0; i < n; i++) cout << " " << ss[i];
		for (i = 0; i < n - 1; i++) {
			custo += D[ ss[i] ][ ss[i+1] ];
		}
		custo += D[ ss[n-1] ][ ss[0] ]; 

		if (custo<custo_opt) {
			custo_opt=custo;
			for (i = 0; i < n; i++) 
				sopt[i]=ss[i];
		}
	}

	int vizinho(int k, int p) // Determina o Vizinho mais próximo do nó p
	{
		p -= 1;
		int i, j; 
		double min = INT_MAX;

		for(i = 0; i < n; i++) {
			if ((k!=i) && (p!=i)) {
				if ((D[k][i]<min) && (set[i]!=1)) {
					min= D[k][i];
					j=i;
				}
			}
		}
		return j;
	}

	void HT()
	{
	   int i, j, c1, c2, i1, i2, cont;

	   int *s1 = new int[n];
	   //AQUI ELE USA A HEURISTICA DO VIZINHO MAIS PROXIMO PARA TODA CIDADE COMO PONTO INICIAL E PEGA O MELHOR
	   for(i=0; i<n; i++) {
		  for(j=0; j<n; j++){
			  set[j]=0; s1[j]=0;
		  }
		  c1=i;
		  set[i]=1;
		  s1[0]=i;
		  // verifica quem é o vert mais perto de i
		  c2= vizinho(i, 0);
		  s1[1]=c2;
		  set[c2]=1;
		  cont=1;
		  // Faz o mesmo para as cabeças c1 e c2 até cont=n
		  while (cont < n -1) {
			 cont++;
			 i1=vizinho(c1, 0);
			 i2=vizinho(c2, 0);
			 if (D[c1][i1]<D[c2][i2]) {
				c1=i1; 
				s1[cont]=i1; 
				set[i1]=1;
			 } else {
				c2=i2; 
				s1[cont]=i2; 
				set[i2]=1;
			 }
		  } // end while

		  g(s1);
	   } // end i
	   delete s1;
	}

	void Resultado(int z)
	{
		int i;
		//cout << endl << custo_opt << endl;
		if (z==1) {
		  cout << "  Rota  = ";
		  for (i = 0; i < n; i++) {
			cout << " " << sopt[i];
		  }
		  cout << endl;
		}
	}

	void encontrar(int origem, int destino, int ss[]);
	

	void vizinhanca(int s[])
	{
		int i, j, k, l, a, b, c, d;

		i=0; j=1; k=2; l=3;
		while (l<n) {
			a=s[i]; b=s[j]; c=s[k]; d=s[l];
			// Avalia a solucao inicial e as soluçoes perturbadas
			g(s);
			// 2a permutação (ijlk), avalia e se desfaz dela
			s[k]=d; s[l]=c; g(s); s[k]=c; s[l]=d;
			// 3a permutação (ikjl), avalia e se desfaz dela
			s[j]=c; s[k]=b; g(s); s[j]=b; s[k]=c;
			// 4a permutação (iklj), avalia e se desfaz dela
			s[j]=c; s[k]=d; s[l]=b; g(s); s[j]=b; s[k]=c; s[l]=d;
			// 5a permutação (iljk), avalia e se desfaz dela
			s[j]=d; s[k]=b; s[l]=c; g(s); s[j]=b; s[k]=c; s[l]=d;
			// 6a permutação (ilkj), avalia e se desfaz dela
			s[j]=d; s[l]=b; g(s); s[j]=b; s[l]=d;
			// 7a permutação (jikl), avalia e se desfaz dela
			s[i]=b; s[j]=a; g(s); s[i]=a; s[j]=b;
			// 8a permutação (jilk), avalia e se desfaz dela
			s[i]=b; s[j]=a; s[k]=d; s[l]=c; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 9a permutação (jkil), avalia e se desfaz dela
			s[i]=b; s[j]=c; s[k]=a; g(s); s[i]=a; s[j]=b; s[k]=c;
			// 10a permutação (jkli), avalia e se desfaz dela
			s[i]=b; s[j]=c; s[k]=d; s[l]=a; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 11a permutação (jlik), avalia e se desfaz dela
			s[i]=b; s[j]=d; s[k]=a; s[l]=c; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 12a permutação (jlki), avalia e se desfaz dela
			s[i]=b; s[j]=d; s[l]=a; g(s); s[i]=a; s[j]=b; s[l]=d;
			// 13a permutação (kijl), avalia e se desfaz dela
			s[i]=c; s[j]=a; s[k]=b; g(s); s[i]=a; s[j]=b; s[k]=c;
			// 14a permutação (kilj), avalia e se desfaz dela
			s[i]=c; s[j]=a; s[k]=d; s[l]=b; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 15a permutação (kjil), avalia e se desfaz dela
			s[i]=c; s[k]=a; g(s); s[i]=a; s[k]=c;
			// 16a permutação (kjli), avalia e se desfaz dela
			s[i]=c; s[k]=d; s[l]=a; g(s); s[i]=a; s[k]=c; s[l]=d;
			// 17a permutação (klij), avalia e se desfaz dela
			s[i]=c; s[j]=d; s[k]=a; s[l]=b; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 18a permutação (klji), avalia e se desfaz dela
			s[i]=c; s[j]=d; s[k]=b; s[l]=a; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 19a permutação (lijk), avalia e se desfaz dela
			s[i]=d; s[j]=a; s[k]=b; s[l]=c; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 20a permutação (likj), avalia e se desfaz dela
			s[i]=d; s[j]=a; s[l]=b; g(s); s[i]=a; s[j]=b; s[l]=d;
			// 21a permutação (ljik), avalia e se desfaz dela
			s[i]=d; s[k]=a; s[l]=c; g(s); s[i]=a; s[k]=c; s[l]=d;
			// 22a permutação (ljki), avalia e se desfaz dela
			s[i]=d; s[l]=a; g(s); s[i]=a; s[l]=d;
			// 23a permutação (lkij), avalia e se desfaz dela
			s[i]=d; s[j]=c; s[k]=a; s[l]=b; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;
			// 24a permutação (lkji), avalia e se desfaz dela
			s[i]=d; s[j]=c; s[k]=b; s[l]=a; g(s); s[i]=a; s[j]=b; s[k]=c; s[l]=d;

			i=i+4; j=j+4; k=k+4; l=l+4;
		} // end while l<=n
	}


	bool g2(int ss[], double &z) {
		int i;
		double custo = 0;
		//for (i = 0; i < n; i++) cout << " " << ss[i];
		for (i = 0; i < n - 1; i++) {
			custo += D[ss[i]][ss[i + 1]];
		}
		custo += D[ss[n - 1]][ss[0]];

		if (custo<z) {
			z = custo;
			for (i = 0; i < n; i++)
				sopt[i] = ss[i];
			return true;
		}
		return false;
	}

	int* vizinhancaz(int s[], double &z)
	{
		int i, j, k, l, a, b, c, d;
		int *sviz = new int[n];
		for (int y = 0; y < n; y++) sviz[y] = s[y];
		i = 0; j = 1; k = 2; l = 3;
		while (l<n) {
			a = s[i]; b = s[j]; c = s[k]; d = s[l];
			// Avalia a solucao inicial e as soluçoes perturbadas
			if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}
			// 2a permutação (ijlk), avalia e se desfaz dela
			s[k] = d; s[l] = c; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[k] = c; s[l] = d;
			// 3a permutação (ikjl), avalia e se desfaz dela
			s[j] = c; s[k] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[j] = b; s[k] = c;
			// 4a permutação (iklj), avalia e se desfaz dela
			s[j] = c; s[k] = d; s[l] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[j] = b; s[k] = c; s[l] = d;
			// 5a permutação (iljk), avalia e se desfaz dela
			s[j] = d; s[k] = b; s[l] = c; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[j] = b; s[k] = c; s[l] = d;
			// 6a permutação (ilkj), avalia e se desfaz dela
			s[j] = d; s[l] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[j] = b; s[l] = d;
			// 7a permutação (jikl), avalia e se desfaz dela
			s[i] = b; s[j] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b;
			// 8a permutação (jilk), avalia e se desfaz dela
			s[i] = b; s[j] = a; s[k] = d; s[l] = c; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 9a permutação (jkil), avalia e se desfaz dela
			s[i] = b; s[j] = c; s[k] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c;
			// 10a permutação (jkli), avalia e se desfaz dela
			s[i] = b; s[j] = c; s[k] = d; s[l] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 11a permutação (jlik), avalia e se desfaz dela
			s[i] = b; s[j] = d; s[k] = a; s[l] = c; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 12a permutação (jlki), avalia e se desfaz dela
			s[i] = b; s[j] = d; s[l] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[l] = d;
			// 13a permutação (kijl), avalia e se desfaz dela
			s[i] = c; s[j] = a; s[k] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c;
			// 14a permutação (kilj), avalia e se desfaz dela
			s[i] = c; s[j] = a; s[k] = d; s[l] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 15a permutação (kjil), avalia e se desfaz dela
			s[i] = c; s[k] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[k] = c;
			// 16a permutação (kjli), avalia e se desfaz dela
			s[i] = c; s[k] = d; s[l] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[k] = c; s[l] = d;
			// 17a permutação (klij), avalia e se desfaz dela
			s[i] = c; s[j] = d; s[k] = a; s[l] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 18a permutação (klji), avalia e se desfaz dela
			s[i] = c; s[j] = d; s[k] = b; s[l] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 19a permutação (lijk), avalia e se desfaz dela
			s[i] = d; s[j] = a; s[k] = b; s[l] = c; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 20a permutação (likj), avalia e se desfaz dela
			s[i] = d; s[j] = a; s[l] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[l] = d;
			// 21a permutação (ljik), avalia e se desfaz dela
			s[i] = d; s[k] = a; s[l] = c; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[k] = c; s[l] = d;
			// 22a permutação (ljki), avalia e se desfaz dela
			s[i] = d; s[l] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[l] = d;
			// 23a permutação (lkij), avalia e se desfaz dela
			s[i] = d; s[j] = c; s[k] = a; s[l] = b; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
			// 24a permutação (lkji), avalia e se desfaz dela
			s[i] = d; s[j] = c; s[k] = b; s[l] = a; if (g2(s, z)) {for (int y = 0; y < n; y++) sviz[y] = s[y];}(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;

			i = i + 4; j = j + 4; k = k + 4; l = l + 4;
		} // end while l<=n
		return sviz;
	}
};

double TSP_instance::f(int ss[])
{
	int i;
	double custo = 0;
	//for (i = 0; i < n; i++) cout << " " << ss[i];
	for (i = 0; i < n - 1; i++) {
		custo += D[ss[i]][ss[i + 1]];
	}
	return custo += D[ss[n - 1]][ss[0]];
}

int* TSP_instance::sol_aleatoria(int s[]) {
	//escolhe duas cidades aleatorias
	int cidade1 = rand() % (n - 1);
	int cidade2 = rand() % (n - 1);
	while (cidade1 == cidade2)
		cidade2 = rand() % (n - 1);

	int *solucao = new int[n];
	for (int i = 0; i < n; i++)
		solucao[i] = s[i];
	//troca elas de posicao
	solucao[cidade1] = s[cidade2];
	solucao[cidade2] = s[cidade1];
	return solucao;
}


void TSP_instance::SimulatedA() {
	double T, alpha, delta, solu_atual, solu_viz;
	int Kint, Kext, p, q;
	int *sol = new int[n];
	int *aux = new int[n];
	double fviz;
	alpha = .7;
	Kint = n;
	Kext = n;
	//iniciar solucao com a do vizinho mais proximo

	HT();
	for (int i = 0; i < n; i++)
		sol[i] = sopt[i];

	T = 100;
	srand(time(NULL));

	p = 1; //iteracoes externas
	double media1, media2;
	do
	{
		media1 = custo_opt;
		cout << "Iter " << p << " Temp" << T << " Cst "<< custo_opt << endl;
		q = 1; //iteracoes internas
		do
		{
			solu_atual = f(sol);
			aux = sol_aleatoria(sol);
			solu_viz = f(aux);
			delta = solu_viz - solu_atual;
			if ((delta < 0) || (rand() < exp(-delta/T))) {
				solu_atual = solu_viz;
				for (int j = 0; j < n; j++)
					sol[j] = aux[j];
				if (solu_atual < custo_opt) {
					for (int j = 0; j < n; j++)
						sopt[j] = sol[j];
					custo_opt = solu_atual;
				}
			}
			q++;
		} while (q < Kint);
		T = alpha*T;
		p++;
		opt4(sopt);
		media2 = custo_opt;
		if (media1 == media2)
			break;
		for (int j = 0; j < n; j++) sol[j] = sopt[n - j - 1]; //perturba
	} while (T > 1);
	
	delete[]sol;
}

void TSP_instance::opt4(int ss[])
{
	int i, j, k, l, a, b, c, d;
	int *s = new int[n];

	for (i = 0; i < n; i++) s[i] = ss[i];
	for (i = 0; i < (n - 3); i++) {
		for (j = i + 1; j < (n - 2); j++) {
			for (k = j + 1; k < (n - 1); k++) {
				for (l = k + 1; l < n; l++) {
					a = s[i]; b = s[j]; c = s[k]; d = s[l];

					// Avalia a 1a permutação ijkl
					g(s);
					// 2a permutação (ijlk), avalia e se desfaz dela
					s[k] = d; s[l] = c; g(s); s[k] = c; s[l] = d;
					// 3a permutação (ikjl), avalia e se desfaz dela
					s[j] = c; s[k] = b; g(s); s[j] = b; s[k] = c;
					// 4a permutação (iklj), avalia e se desfaz dela
					s[j] = c; s[k] = d; s[l] = b; g(s); s[j] = b; s[k] = c; s[l] = d;
					// 5a permutação (iljk), avalia e se desfaz dela
					s[j] = d; s[k] = b; s[l] = c; g(s); s[j] = b; s[k] = c; s[l] = d;
					// 6a permutação (ilkj), avalia e se desfaz dela
					s[j] = d; s[l] = b; g(s); s[j] = b; s[l] = d;
					// 7a permutação (jikl), avalia e se desfaz dela
					s[i] = b; s[j] = a; g(s); s[i] = a; s[j] = b;
					// 8a permutação (jilk), avalia e se desfaz dela
					s[i] = b; s[j] = a; s[k] = d; s[l] = c; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 9a permutação (jkil), avalia e se desfaz dela
					s[i] = b; s[j] = c; s[k] = a; g(s); s[i] = a; s[j] = b; s[k] = c;
					// 10a permutação (jkli), avalia e se desfaz dela
					s[i] = b; s[j] = c; s[k] = d; s[l] = a; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 11a permutação (jlik), avalia e se desfaz dela
					s[i] = b; s[j] = d; s[k] = a; s[l] = c; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 12a permutação (jlki), avalia e se desfaz dela
					s[i] = b; s[j] = d; s[l] = a; g(s); s[i] = a; s[j] = b; s[l] = d;
					// 13a permutação (kijl), avalia e se desfaz dela
					s[i] = c; s[j] = a; s[k] = b; g(s); s[i] = a; s[j] = b; s[k] = c;
					// 14a permutação (kilj), avalia e se desfaz dela
					s[i] = c; s[j] = a; s[k] = d; s[l] = b; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 15a permutação (kjil), avalia e se desfaz dela
					s[i] = c; s[k] = a; g(s); s[i] = a; s[k] = c;
					// 16a permutação (kjli), avalia e se desfaz dela
					s[i] = c; s[k] = d; s[l] = a; g(s); s[i] = a; s[k] = c; s[l] = d;
					// 17a permutação (klij), avalia e se desfaz dela
					s[i] = c; s[j] = d; s[k] = a; s[l] = b; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 18a permutação (klji), avalia e se desfaz dela
					s[i] = c; s[j] = d; s[k] = b; s[l] = a; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 19a permutação (lijk), avalia e se desfaz dela
					s[i] = d; s[j] = a; s[k] = b; s[l] = c; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 20a permutação (likj), avalia e se desfaz dela
					s[i] = d; s[j] = a; s[l] = b; g(s); s[i] = a; s[j] = b; s[l] = d;
					// 21a permutação (ljik), avalia e se desfaz dela
					s[i] = d; s[k] = a; s[l] = c; g(s); s[i] = a; s[k] = c; s[l] = d;
					// 22a permutação (ljki), avalia e se desfaz dela
					s[i] = d; s[l] = a; g(s); s[i] = a; s[l] = d;
					// 23a permutação (lkij), avalia e se desfaz dela
					s[i] = d; s[j] = c; s[k] = a; s[l] = b; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;
					// 24a permutação (lkji), avalia e se desfaz dela
					s[i] = d; s[j] = c; s[k] = b; s[l] = a; g(s); s[i] = a; s[j] = b; s[k] = c; s[l] = d;

					l = l + salto;
				} // end for l
				k = k + salto;
			} // end for k
			j = j + salto;
		} // end for j
		i = i + salto;
	} // end for i
}