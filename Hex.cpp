#include<iostream>
#include<string>
#include<vector>
#include<time.h>
#include<dos.h>
#include<Windows.h>
//#include<synchapi.h>
#include<set>


using namespace std;


#define Assert(X,Y) if(!(X)){printf(Y);while(1);}

#define IDX(r,c) (((r)*gnBS)+(c))
#define GET_COL(idx) ((idx)%gnBS)
#define GET_ROW(idx) ((idx)/gnBS)

int gnBS;

class BoardDisplay
{
	string* BoardChars;
	int BoardSize;

public:
	BoardDisplay(int n) :BoardSize(n)
	{
		BoardChars = new string[n];
		string Pattern;
		for (int i = 0; i < n; ++i)
			Pattern.push_back('.');
		for (int i = 0; i < n; ++i)
			BoardChars[i] = Pattern;  //.resize(n);
		//BoardChars[1][1] = 'X';
		//BoardChars[2][2] = 'O';
	}
	string* getBoardChars(){ return BoardChars; }

	void printDots(string str)
	{
		int ss = str.size();
		for (int i = 0; i < ss; ++i)
		{
			cout << str[i];
			if (i == ss - 1)continue;
			cout << " - ";
		}
	}

	void printSlashes(int n)
	{
		while (n--)
		{
			cout << "\\";
			if (n == 0)continue;
			cout << " / ";
		}
	}


	void Display_PrintBoard()
	{
		int rows = 2 * BoardSize - 1;
		cout << "  ";
		for (int i = 0; i < BoardSize; ++i)
		{
			cout<< i << "   ";
		}
		cout << endl;
		for (int i = 0; i < rows; ++i)
		{
			int space = i;
			while (space--)cout << " ";
			if (i % 2 == 0)
			{
				cout <<"   "<<"\b\b"<<i / 2<<" ";
				printDots(BoardChars[i / 2]);
			}
			else
			{
				cout << "   ";
				printSlashes(BoardSize);
			}
			
			cout << endl;
		}
	}

	bool Display_MakeMove(int r,int c,char ch)
	{
		Assert((r >= 0 && r < BoardSize), "Wrong Index to make move");
		Assert((c >= 0 && c < BoardSize), "Wrong Index to make move");
		Assert((ch == 'X' || ch == 'O'), "Wrong Character");
		if (BoardChars[r][c] == '.')
		{
			BoardChars[r][c] = ch;
			return true;
		}
		return false;
	}
	void reset()
	{
		for (int i = 0; i < BoardSize; ++i)
			for (int j = 0; j < BoardSize; ++j)
				BoardChars[i][j] = '.';
	}
};

class Graph{
	int V,E;
	vector<int> *Adj;

public:
	Graph(int v, int e=0) :V(v),E(e){
		Adj = new vector<int>[V];
	}
	void addEdge(int a, int b)
	{
		Adj[a].push_back(b);
		E++;
	}
	void addUndirectedEdge(int a, int b)
	{
		Assert(a < V && b < V, "Wrong Node Number");
		Adj[a].push_back(b);
		Adj[b].push_back(a);
		E += 2;
	}
	//int PrimMST();
	vector< int> * getAdjList(){ return Adj; }
	int getEdges(){ return E; }
};



class Dijstra
{
	Graph* g;
	multiset<pair<int, int> >pq;
	bool *vis;
	int *d;
	int * prev;
	int Nodes;
public:
	Dijstra(int v, Graph* G) :g(G), Nodes(v)
	{
		vis = new bool[v]();
		d = new int[v]();
		prev = new int[v]();
		for (int i = 0; i < v; ++i)
		{
			prev[i] = -1;
			d[i] = 1e9;
			vis[i] = false;
		}
	}
	void Init()
	{
		for (int i = 0; i < Nodes; ++i)
		{
			d[i] = 1e9;
			prev[i] = -1;
			vis[i] = false;
		}
	}

	void calShortestPath(int s, char ch,const char* chArr)
	{
		vector<int > *Adj = g->getAdjList();
		d[s] = 0;
		pq.insert(make_pair(0, s));
		prev[s] = s;
		while (!pq.empty())
		{
			pair<int, int> p = *pq.begin();
			pq.erase(pq.begin());
			int cn = p.second;
			int cw = p.first;
			if (vis[cn])continue;
			vis[cn] = true;

			for (int i = 0; i < Adj[cn].size(); ++i)
			{
				int pp = Adj[cn][i];
				int node = pp;
				if (chArr[node] != ch)
					continue;
				int wgt = 1;
				if (d[cn] + wgt < d[node])
				{
					d[node] = d[cn] + wgt;
					pq.insert(make_pair(d[node], node));
					prev[node] = cn;
				}
			}

		}
	}


	int * getDistanceArray(){ return d; }
	int * getPathArray(){ return prev; }
	bool* getVis(){ return vis; }
};

#define DOUBLELOOP for(int i=0;i<BoardSize;++i) \
					for (int j = 0; j < BoardSize;++j)


int gDbg = 0;
class MonteCarloSim{
	int * EFA;
	int *curEFA;
	int ep;
	int nN;
	int cep;
public:
	MonteCarloSim(char * BOrig, int V)
	{
		nN = V;
		EFA = new int[V];
		curEFA = new int[V];
		ep = V;
		cep = V;
		for (int i = 0; i < ep; ++i)
		{
			
			EFA[i] = i;
			if (BOrig[i] != '.')
			{
				while (BOrig[ep - 1] != '.')
					ep--;
				if (i >= ep)continue;

				swap(EFA[i], EFA[--ep]);
				EFA[i] = ep;

			}
		}
	}

	void newMCRun(char* BOrig)
	{
		cep = nN;
	//	for (int i = 0; i < cep; ++i)
		//	curEFA[i] = i;
		for (int i = 0; i < cep; ++i)
		{
			curEFA[i] = i;
			if (BOrig[i] != '.')
			{
				while (BOrig[cep - 1] != '.')
				{
					cep--;
					if (i >= cep)break;
				}
				if (i >= cep)continue;
				swap(curEFA[i], curEFA[--cep]);

				curEFA[i] = cep;
			}
		}
	}

	int getRandomIdx()
	{
		if (ep == 0)return -1;
		int rIdx = rand() % ep;
		int idx = EFA[rIdx];
		swap(EFA[rIdx], EFA[--ep]);
		return idx;
	}
	int getRandomIdxforCurRun()
	{
		if (cep == 0)return -1;
		int rIdx = rand() % cep;
		int idx = curEFA[rIdx];
		swap(curEFA[rIdx], curEFA[--cep]);
		return idx;
	}

	void updateEFAwithMove(int idx)
	{
		swap(EFA[idx], EFA[--ep]);
	}

	~MonteCarloSim()
	{
		delete[] EFA;
		delete[] curEFA;
	}

};


class HexBoard{
	BoardDisplay BD;
	int BoardSize;
	Graph g;
	char * chArr;
	//char * EmptyFirst;
	int OccupiedNodes;
	bool *vis;
	Dijstra *Dj;
	int nNodes;
public:
	HexBoard(int n) :BoardSize(n), BD(n), g(n*n), OccupiedNodes(0), nNodes(n*n)
	{
		gnBS = n;
		chArr = new char[n*n];
		//EmptyFirst = new char[n*n];
		
		for (int i = 0; i < n*n; ++i)
			chArr[i] = '.';
		LinkNodes();
		vis = new bool[n*n];
		Dj = new Dijstra(BoardSize*BoardSize, &g);

	}

	int printPath(int *p, int idx)
	{
		if (p[idx] == idx) return idx;
		return printPath(p, p[idx]);
		//cout << "->" << p[idx];
	}

	void getLargestPath(char ch,int *src,int *dst)
	{
		Dj->Init();
		bool* DVis = Dj->getVis();
		DOUBLELOOP
		{
			int idx = IDX(i, j);
			if (!DVis[idx] && chArr[idx] == ch)
				Dj->calShortestPath(idx, ch,chArr);
		}
		
		int * DArr = Dj->getDistanceArray();
		int maxDist = INT_MIN;
		for (int i = 0; i < BoardSize*BoardSize; ++i)
		{
			if (DArr[i] < 1e9 && DArr[i] > maxDist)
			{
				maxDist = DArr[i];
				*dst = i;
			}
		}
		*src = printPath(Dj->getPathArray(), *dst);

	}


	bool stopLongestPath()
	{
		int s, d;
		getLargestPath('X', &s, &d);
		if (GET_COL(s)>GET_COL(d))
			swap(s, d); //s always on left

		bool bMoveMade = false;
		//try make move at left side
		if (GET_COL(s) > 0)
		{
			vector<int>* NodePtr = g.getAdjList() + s;

			for (int i = 0; i < NodePtr->size(); ++i)
			{
				int it = NodePtr->at(i);
				if (GET_COL(it)<GET_COL(s) && chArr[it] == '.')
				{
					bMoveMade = MakeMove(GET_ROW(it), GET_COL(it), 'O');
					Assert(bMoveMade, "Move not Made");
					break;
				}
			}
		}
		if (!bMoveMade && GET_COL(d) < BoardSize - 1)
		{
			vector<int>* NodePtr = g.getAdjList() + d;

			for (int i = 0; i < NodePtr->size(); ++i)
			{
				int it = NodePtr->at(i);
				if (GET_COL(it)>GET_COL(d) && chArr[it] == '.')
				{
					bMoveMade = MakeMove(GET_ROW(it), GET_COL(it), 'O');
					Assert(bMoveMade, "Move not Made");
					break;
				}
			}
		}
		if (!bMoveMade)
			RandomMove('O');
		return bMoveMade;
	}

	void LinkHorizontalNodes(int r,int c)
	{
		if (c >= BoardSize - 1)return;
		g.addUndirectedEdge(IDX(r, c), IDX(r, c + 1));
	}

	void LinkVerticalNodes(int r, int c)
	{
		if (r >= BoardSize - 1)return;
		g.addUndirectedEdge(IDX(r, c), IDX(r+1, c));
	}

	void LinkBackwayNodes(int r, int c)
	{
		if (r >= BoardSize-1 || c <= 0 )return;
		g.addUndirectedEdge(IDX(r, c), IDX(r + 1, c - 1));
	}

	int CheckEdgesCount(int n)
	{
		int thEdges = (2 * ((n - 1)*n) + (n - 1)*(n - 1))*2;
		int ActualEdges = g.getEdges();
		Assert(thEdges == ActualEdges, "Wrong Number of Edges");
		return ActualEdges;
	}

	void LinkNodes()
	{

		for (int i = 0; i < BoardSize; ++i)
		{
			for (int j = 0; j < BoardSize; ++j)
			{
				LinkHorizontalNodes(i,j);
				LinkVerticalNodes(i,j);
				LinkBackwayNodes(i,j);
			}
		}
		CheckEdgesCount(BoardSize);
		
	}
	void PrintBoard()
	{
		BD.Display_PrintBoard();
	}
	bool MakeMove(int r,int c,char ch)
	{
		Assert(OccupiedNodes < BoardSize*BoardSize, "Board Full");
		if (r >= BoardSize || c >= BoardSize)return false;
		if (r < 0 || c < 0) return false;
		if (ch != 'X' && ch != 'O')return false;
		if (BD.Display_MakeMove(r, c, ch))
		{
			//Assert(chArr[IDX(r, c)])
			chArr[IDX(r, c)] = ch;
			OccupiedNodes++;
			return true;
		}
		return false;
	}

	void findEmptyNode(int *r, int*c)
	{
		for (int i = 0; i < BoardSize; ++i)
		{
			for (int j = 0; j < BoardSize; ++j)
			{
				if (chArr[IDX(i, j)] == '.')
				{
					*r = i;
					*c = j;
					return;
				}
			}
		}
		Assert(false, "No Empty Nodes found");
	}

	void RandomMove(char ch)
	{
		int retryCnt = BoardSize*BoardSize;
		
		while (!MakeMove(rand() % BoardSize, rand() % BoardSize, ch) && retryCnt>0)
			retryCnt--;
		if (retryCnt == 0)
		{
			int r, c;
			findEmptyNode(&r,&c);
			Assert(MakeMove(r, c,ch), "Invalid Empty Node");
		}
	}

	bool isBoardFull()
	{
		return (OccupiedNodes == BoardSize*BoardSize);
	}
	
	bool DFS(int r, int c, char ch)
	{
		int idx = IDX(r, c);
		vis[idx] = true;
		vector<int>* NodePtr = g.getAdjList() + idx;
		for (int i = 0; i < NodePtr->size(); ++i)
		{
			int it = NodePtr->at(i); 
			if (!vis[it] && chArr[it] == ch)
			{
				if (ch == 'X' && GET_COL(it) == BoardSize - 1)
					return true;
				if (ch == 'O' && GET_ROW(it) == BoardSize - 1)
					return true;
				if (DFS(GET_ROW(it), GET_COL(it), ch))
					return true;
			}
		}
		return false;
	}

	bool CheckP1Winner(char ch)
	{
		memset(vis, false, BoardSize*BoardSize);
		//bool winP1 = false;
		for (int i = 0; i < BoardSize; ++i)
		{
			int idx = IDX(i, 0);
			if (!vis[idx] && chArr[idx] == ch)
			{
				if (DFS(i, 0, ch))
					return true;
			}
		}

		return false;
	}

	bool CheckP2Winner(char ch)
	{
		memset(vis, false, BoardSize*BoardSize);
		//bool winP1 = false;
		for (int i = 0; i < BoardSize; ++i)
		{
			if (!vis[IDX(0, i)] && chArr[IDX(0, i)] == ch)
			{
				if (DFS(0, i, ch))
					return true;
			}
		}

		return false;
	}

	void P2_blockP1Move()
	{
		for (int j = BoardSize - 1; j >= 0; j--)
		{
			for (int i = 0; i < BoardSize; ++i)
			{
				int idx = IDX(i, j);
				vector<int>* NodePtr = g.getAdjList() + idx;
				for (int k = 0; k < NodePtr->size(); ++k)
				{
					int it = NodePtr->at(k);
					if (GET_COL(it) < j && chArr[it]=='X')
					{
						if(MakeMove(i, j, 'O'))
						return;
					}
				}
			}
		}
		int r, c;
		findEmptyNode(&r, &c);
		Assert(MakeMove(r, c, 'O'), "Invalid Empty Node");
	}
	const int MAX_SIM = 1000;
	const int MAX_P = 2;

	void rebuildFromBkupArr(char* Bkup)
	{
		OccupiedNodes = 0;
		BD.reset();
		memset(chArr, '.', nNodes);
		for (int j = 0; j < nNodes; ++j)
		{
			if (Bkup[j] != '.')
				Assert(MakeMove(GET_ROW(j), GET_COL(j), Bkup[j]), "Wrong Move");
		}
	}
	int MakeMonteCarloMove(char ch)
	{
		//make copy of chArr
		char * Bkup = new char[nNodes];
		memcpy(Bkup, chArr, nNodes);
		int BkupON = OccupiedNodes;

		int * WinArr = new int[nNodes]();

		//create MC Object
		//create Empty First Arr from chArr
		MonteCarloSim MCSim(chArr, nNodes);

		//for all empty Nodes
		int RandIdx = MCSim.getRandomIdx();
			//get one random point from Empty First Arr
		while (RandIdx >= 0)
		{
		
			//update Empty First Arr
			//Run 1000 
			for (int i = 0; i < MAX_SIM;++i)
			{
				Assert(MakeMove(GET_ROW(RandIdx), GET_COL(RandIdx), ch), "Wrong Move");
			//	MCSim.updateEFAwithMove(RandIdx);
				//init MC Sim
				MCSim.newMCRun(chArr);
				int curRandIdx = MCSim.getRandomIdxforCurRun();
				char turns[2] = { 'X', 'O' };
				int tIdx = 0;
				//while BoardNotFull
				while (curRandIdx >= 0)
				{
					//Get Random for Cur
					//Make Move
					Assert(MakeMove(GET_ROW(curRandIdx), GET_COL(curRandIdx), turns[tIdx++%2]), "Wrong Move");

					curRandIdx = MCSim.getRandomIdxforCurRun();
				}

				Assert(OccupiedNodes == nNodes, "Board Not Filled");

				//Check Winner and Update Win Arr
				if (CheckP2Winner('O'))
					WinArr[RandIdx]++;
				else
					Assert(CheckP1Winner('X'), "No Winner");

				//cout << "---------------------Board After Simul------------"<<endl;
				//PrintBoard();
				//rebuild from BKUP chArr
				rebuildFromBkupArr(Bkup);
				
			}
			RandIdx = MCSim.getRandomIdx();
		}
		//Choose optimum
		int mi = 0;
		for (int k = 1; k < nNodes; ++k)
		{
			if (WinArr[k]>WinArr[mi])mi = k;
		}
		//Rebuild from chArr
		rebuildFromBkupArr(Bkup);
		//Make Move
		Assert(MakeMove(GET_ROW(mi), GET_COL(mi), ch), "Wrong Move");
		return mi;
	}
};

#define P1_PLAYER 1
#define AUTO_DEBUG 0
class GamePlay
{
	HexBoard HB;
	char P1, P2;
public:
	GamePlay(int n) :HB(n), P1('X'), P2('O')
	{
		HB.PrintBoard();
	}

	void GameLoop()
	{
		do{
#if(P1_PLAYER)

			int r, c;
			bool validMove = false;
			do{
				cout << "Enter Your Move in form of row col : ";
#if AUTO_DEBUG
				r = 2;c=2;
#else
				cin >> r >> c;
#endif
				validMove = HB.MakeMove(r, c, P1);
				if (!validMove)
					cout << "Wrong Move: Enter row col again" << endl;
				else
				{
					cout << "Your Move" << endl << endl;
					HB.PrintBoard();
					//delay(1000);
					Sleep(1000);
				}
			} while (!validMove);
#else
			HB.RandomMove('X');
			cout << "Your Move" << endl << endl;
			HB.PrintBoard();
			//delay(1000);
			Sleep(500);
#endif
			if (HB.CheckP1Winner('X'))
			{
				cout << "---------Congrats Player Won---------";
				break;
			}
			if (HB.isBoardFull())
			{
				cout << "------------Board Full------------";
				break;
			}

			//if (rand() % 3)
				//HB.stopLongestPath();
			int cm =	HB.MakeMonteCarloMove('O');
				//HB.P2_blockP1Move();
			//else
				//HB.RandomMove(P2);
			cout << "Computer's Move r" <<GET_ROW(cm)<<" c"<<GET_COL(cm)<< endl;

			HB.PrintBoard();
			if (HB.CheckP2Winner('O'))
			{
				cout << "---------Sorry Computer Won---------";
				break;
			}
		} while (1);
	}
};


int main()
{
	srand(time(NULL));
	//srand(12345);
	int n; cout << "Enter Board Size: ";
#if AUTO_DEBUG
	n=5;
#else
	cin >> n;
#endif
	cout << "Initial Board :" << endl << endl;
	GamePlay GP(n);
	GP.GameLoop();
	//HexBoard HB(n);
	//HB.PrintBoard();
	while (!(cin >> n));
}