#include <iostream>
#include <time.h>
#include "heuristix.h"
#include <algorithm>
using namespace std;
#define IN 99
void timeused(double *time);

void TSP_instance::HeuristicaRUIM(){

	int i, j;

	int *s1 = new int[n];
/**
	Ideia: pegar vertices com menores desvios padroes
	
	int minV;
	double *DesvP = new double[n];
	//double *distg = new double[n];
	double *E = new double[n];
	for(i=0; i<n; i++){
		set[i] = 0;
		E[i] = 0;
		for(j=0; j < n; j++){
			if(i != j)
				E[i] += D[i][j];
		}
		E[i] = E[i]/(n-1);
		double soma = 0;
		for(j=0; j < n; j++){
			soma += pow(D[i][j] - E[i],2)/(n-1);
		}
		DesvP[i] = sqrt(soma);
	}
	cout << endl << endl;
	int imin;
	double iminV;
	for (i = 0; i < n; i++){
		iminV = INT_MAX;
		for (j = 0; j < n; j++){
			if (DesvP[j] < iminV && set[j] == 0){
				imin = j;
				iminV = DesvP[j];
			}
		}
		set[imin] = 1;
		s1[i] = imin;
		//cout << imin << " ";
	}
	g(s1);
	vizinhanca(s1);

	int troca;
	int *s2 = new int[n];
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			for (int w = 0; w < n; w++){
				s2[w] = s1[w];
			}

			if (i != j){
				s2[i] = s1[j];
				s2[j] = s1[i];
				g(s2);
			}
		}
		
	}
	*/
	for(i = 0; i < n; i++)
	{
		encontrar(0,n-1,s1);
		getchar();
	}
	delete s1;
}

void TSP_instance::encontrar(int origem, int destino, int ss[]) {

	// Vetor que guarda a menor distancia conhecida desde a origem ate'
	// cada vertice. Se nenhuma distancia ainda e' conhecida, dizemos
	// que a menor distancia conhecida ate o vertice e' igual a infinito.
	int *melhorDist = new int[n]; // nao esquecer do 'delete'
	for (int i=0; i<n; i++)
		melhorDist[i] = -1; // -1 significa "infinito"

	// Vetor que informa se ja conhecemos a menor distancia desde o
	// vertice de origem ate' um dado vertice.
	bool *concluido = new bool[n]; // nao esquecer do 'delete'
	for (int i=0; i<n; i++)
		concluido[i] = false; // false significa "ainda nao concluido"

	// Informa o predecessor de cada vertice no caminho desde a origem
	int *pred = new int[n]; // nao esquecer do 'delete'
	for (int i=0; i<n; i++)
		pred[i] = -1; // -1 significa "sem antecessor"

	// A busca inicia-se a partir do vertice "origem"
	melhorDist[origem] = 0; // Distancia de "origem" para si mesmo: zero.

	while (1) {

		// Encontrar vertice nao-concluido com menor distancia
		int menor = -1; // Esta variavel armazena o indice de tal vertice
		for (int i=0; i<n; i++) {
			if ((concluido[i] == false) && (melhorDist[i] >= 0)) {
				if (menor == -1)
					menor = i;
				else
					if (melhorDist[i] < melhorDist[menor])
						menor = i;
			}
		}

		// Se "menor" e' o vertice de destino, a tarefa foi concluida!
		if (menor == destino) {
			cout << "Distancia minima de " << origem << " para ";
			cout << destino << " e': " << melhorDist[destino] << endl;
			break;
		}

		// Se "menor" e' -1, nao ha vertice com valor finito de "melhorDist".
		// A conclusao e' que nao temos como chegar ao vertice "destino".
		if (menor == -1) {
			cout << "Nao consigo atingir o vertice de destino." << endl;
			return;
		}

		// Marcar o vertice "menor" como concluido
		concluido[menor] = true;

		// Considerar a atualizacao das melhores distancias conhecidas ate'
		// os vertices que sao diretamente alcancaveis a partir de "menor".
		for (int i=0; i<n; i++) {
			if ((D[menor][i] > 0) && (concluido[i] == false)) {
				// Se nao conhecemos ainda nenhum caminho ate o vertice "i"
				if (melhorDist[i] < 0) {
					melhorDist[i] = melhorDist[menor] + D[menor][i];
					pred[i] = menor;
				}
				else // Ou se o caminho existente e' pior que o descoberto
					if (melhorDist[i] > melhorDist[menor] + D[menor][i]) {
						melhorDist[i] = melhorDist[menor] + D[menor][i];
						pred[i] = menor;
					}
			}
		}

	}// fim do while

	// Imprimir caminho (invertido)
	int k = destino;
	cout << k;
	while (pred[k] != -1) {
		cout << " <- " << pred[k];
		k = pred[k];
	}
	cout << endl;

	// Nao esquecer: desalocar "melhorDist" e "concluido".
	delete melhorDist;
	delete concluido;
	delete pred;
}


int main(int argc, char *argv[])
{
	if(argc < 3){
		cerr << "	Erro: Sem argumentos de entrada suficientes" << endl;
		system("pause");
		exit(0);
	}
	double time;
	ofstream resultados;
	resultados.open("resultados.txt", fstream::app);
	//ATENCAO -- CONVERTER ARGUMENTOS DE CHAR PARA INT USAR FUNÇÃO ''ATOI''
	int tipoM = atoi(argv[2]);

    TSP_instance Problema(argv[1], tipoM);
	/**
    timeused(NULL);
    Problema.HT();
    timeused(&time);
	cout << endl << endl << "A solucao de HHVP gastou " << time << " seg. e foi : " << Problema.custo_opt;
    Problema.Resultado(1);
	*/
	/**
	timeused(NULL);
    Problema.HeuristicaRUIM();
    timeused(&time);
	cout << endl << endl << "A solucao de Heuristica RUIM gastou " << time << " seg. e foi : " << Problema.custo_opt;
    Problema.Resultado(1);
	*/
	timeused(NULL);
	Problema.SimulatedA();
	timeused(&time);
	cout << endl << endl << "A solucao de SA gastou " << time << " seg. e foi : " << Problema.custo_opt;
	Problema.Resultado(1);
	resultados << argv[1] << " " << Problema.custo_opt << " " << time << endl;
	
    resultados.close();
	



    //system("pause");
    return 0;
}

void timeused(double *time)
{
  static double tstart, tend, tprev;

  if (time == NULL) {
    clock(); /* one extra call to initialize clock */
    tstart = tprev = clock();
  } else {
    tend = clock();
    if (tend < tprev) tstart -= ULONG_MAX; /* wraparound occured */
    tprev = tend;
    *time = (tend-tstart) / CLOCKS_PER_SEC; /* convert to seconds */
  }
}

