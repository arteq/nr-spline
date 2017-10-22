/*  PROGRAM: spline
 *  WERSJA: 0.2
 *  Autor: ArteQ: <arteq(at)arteq(dot)org>
 *  OSTATNIA MODYFIKACJA: 2011/11/26 (sob) 16:22:54
 *  KEYWORDS: spline, metody numeryczne, interpolacja, funkjce sklejane, cubic
 *		spline interpolation, rozkład sygnału na mody empityczne, emd, csi
 */

#include	<iostream>
#include	<fstream>
#include	<math.h>
#include	<vector>
#include	"nrutil.h"

using namespace std;

int ile;	// ilosc punktow do wczytania
int nr_max = 0, nr_min = 0;


/* ====================================================================== */

void wczytaj_dane(char plik[], float *dat_x, float *dat_y)
{
	ifstream fp;
	int i;
	float x,y;

	fp.open(plik, ios::in);
	for (i = 1; i < ile+1; i++)
	{
		fp >> x;
		fp >> y;

		dat_x[i] = x;
		dat_y[i] = y;
	}

	fp.close();
	cout << "Wczytano dane do pamieci..." << endl;
}

/* ====================================================================== */

int rozmiar(char plik[])
{
	ifstream fp;
	int i=0;
	float x,y;

	fp.open(plik, ios::in);
	while ( !fp.eof() )
	{
		fp >> x;
		fp >> y;

		i++;
	}
	fp.close();

	i--;
	cout << "# Plik zawiera: " << i << " punktow danych\n" << endl;
	return i;
}


/* ====================================================================== */

/* x[] - dane wejsciowe - x
 * y[] - dane wejsciowe - y
 * n - liczba danych ktore lykamy
 * yp1, ypN - pierwsze pochodne na poczatku i koncu przedzialu
 * y2[] - odliczone drugie pochodne funkcji
 */
void spline(float *x, float *y, int n, float yp1, float ypn, float *y2)
{
	float p, qn, sig, un, *u;
	u = vectorek(1, n-1);


	cout << "Obliczanie pochodnych..." << endl;
	if (yp1 > 0.99e30)	// splajnik naturalny...
		y2[1] = u[1] = 0.0;
	else								// ...albo i nie :-)
	{
		y2[1] = -0.5;
		u[1] = (3.0 / (x[2] - x[1]) ) * ( (y[2] - y[1]) / (x[2] - x[1]) - yp1 );
	}
	
	int i;
	for (i=2; i<=n-1; i++)
	{
		sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
		u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
	}

	if (ypn > 0.99e30)
		qn = un = 0.0;
	else
	{
		qn = 0.5;
		un = (3.0 / (x[n] - x[n-1]) ) * (ypn - (y[n] - y[n-1]) / (x[n] - x[n-1]) );
	}

	y2[n] = (un - qn * u[n-1]) / (qn * y2[n-1] + 1.0);

	for (int k = n-1; k>=1; k--)
	{
		y2[k] = y2[k] * y2[k+1] + u[k];
	}

	free_vector(u, 1, n-1);
	cout << "Pochodne obliczone..." << endl;
}

/* ====================================================================== */

/* xa[] - funkcja dana X
 * ya[] - funkcja dana Y
 * n - ilosc punktow
 * y2a[] - tablica pochodnych policzona wyzej
 * x - arg dla ktorego zwracany jest...
 * y - igrek :-)
 */
float splint(float *xa, float *ya, float *y2a, int n, float x)
{
	int klo, khi, k;
	float h, b, a;
	float wynik;

	klo = 1;
	khi = n;

	while (khi - klo > 1)
	{
		k = (khi + klo) >> 1;
		if (xa[k] > x) khi = k;
		else klo = k;
	}

	h = xa[khi] - xa[klo];
	if (h == 0) cout << "Wystapil bledny blad :-(" << endl;
	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h;

	wynik = a * ya[klo] + b * ya[khi] + ( (a*a*a - a) * y2a[klo] + (b*b*b - b)*y2a[khi] ) * (h*h) / 6.0;
	return wynik;
}

/* ====================================================================== */

void szukaj_ekstremow(float *daneX, float *daneY, 
		vector<float> *maxX, vector<float> *maxY, vector<float> *minX, vector<float> *minY, int n)
{

	for (int i = 1; i < n; i++)
	{
		if (daneY[i-1] < daneY[i] && daneY[i+1] < daneY[i])
		{
			maxX->push_back( daneX[i] );
			maxY->push_back( daneY[i] );
		}
		if (daneY[i-1] > daneY[i] && daneY[i+1] > daneY[i])
		{
			minX->push_back( daneX[i] );
			minY->push_back( daneY[i] );
		}
	}
}

/* ====================================================================== */

void zapisz_ekstrema(vector<float> &maxX, vector<float> &maxY, vector<float> &minX, vector<float> &minY)
{
	ofstream fp;

	fp.open("max.dat", ios::out);
	for (int i = 0; i < maxX.size(); i++)
		fp << maxX[i] << "\t" << maxY[i] << endl;
	fp.close();

	fp.open("min.dat", ios::out);
	for (int i = 0; i < minX.size(); i++)
		fp << minX[i] << "\t" << minY[i] << endl;
	fp.close();

	cout << "\nMaksimow: " << maxX.size() << "\tMinimow: " << minX.size() << endl;
}

/* ====================================================================== */

void ekstrema_to_tab(vector<float> &maxX, vector<float> &maxY, vector<float> &minX, vector<float> &minY,
		float *dat_max_X, float *dat_max_Y, float *dat_min_X, float *dat_min_Y)
{

	for (int i = 1; i < maxX.size() +1; i++)
	{
		dat_max_X[i] = maxX[i-1];
		dat_max_Y[i] = maxY[i-1];
	}

	for (int i = 1; i < minX.size() +1; i++)
	{
		dat_min_X[i] = minX[i-1];
		dat_min_Y[i] = minY[i-1];
	}
}

/* ====================================================================== */

int main()
{
	float *dat_x, *dat_y, *pochodne, *pochodne_max, *pochodne_min, *spl;
	float *dat_min_X, *dat_min_Y, *dat_max_X, *dat_max_Y;
	vector<float> minX, minY, maxX, maxY;
	char nazwa[80] = "zad2.txt";
	ofstream output;

	ile = rozmiar(nazwa);
	float argument;

	dat_x = new float[ile+1];
	dat_y = new float[ile+1];
	pochodne = new float[ile+1];
	spl = new float[ile+1];

	cout << "Zarezerwowano pamiec..." << endl;

	wczytaj_dane(nazwa, dat_x, dat_y );
	spline(dat_x, dat_y, ile, 0, 0, pochodne);

	output.open("spline.dat", ios::out);
	output << "#X \tY \tpoch \tSPLINE" << endl;

	int interpol = 10000;

	for (int x = 1; x < ile+1; x++)
	{
		spl[x] = splint(dat_x, dat_y, pochodne, ile+1, (1.0*x)/100.0);
		output << dat_x[x] << "\t" << dat_y[x] << "\t" << pochodne[x] << "\t" << spl[x] << endl;
	}

	output.close();

	szukaj_ekstremow(dat_x, dat_y, &maxX, &maxY, &minX, &minY, ile);
	zapisz_ekstrema(maxX, maxY, minX, minY);

	pochodne_max = new float[maxX.size()+1 ];
	pochodne_min = new float[minX.size()+1 ];

	dat_max_X = new float[ maxX.size() + 1];
	dat_max_Y = new float[ maxX.size() + 1];
	dat_min_X = new float[ minX.size() + 1];
	dat_min_Y = new float[ minX.size() + 1];

	ekstrema_to_tab(maxX, maxY, minX, minY, dat_max_X, dat_max_Y, dat_min_X, dat_min_Y);
	spline(dat_max_X, dat_max_Y, maxX.size(), 0, 0, pochodne_max);
	spline(dat_min_X, dat_min_Y, minX.size(), 0, 0, pochodne_min);


	cout << "Mod empiryczny " << endl;
	output.open("mod.dat", ios::out);
	for (int x = 0; x < ile; x++)
	{
		float Y, Z;
		Y = splint(dat_max_X, dat_max_Y, pochodne_max, maxX.size(), (1.0*x)/100.0);
		Z = splint(dat_min_X, dat_min_Y, pochodne_min, minX.size(), (1.0*x)/100.0);
		output << x/100.0 << "\t" << dat_x[x] - 0.5*(Y+Z) << "\n";
	}
	output.close();

	cout << "\n";
	cout << "********************************************************************************\n";
	cout << "* Program zapisal pliki:                                                       *\n";
	cout << "* -> spline.dat      zawierajacy wynik interpolacji CSI                        *\n";
	cout << "* -> max.dat         zawierajacy polozenia i wartosci maximow sygnalu          *\n";
	cout << "* -> min.dat         zawierajacy polozenia i wartosci minimow sygnalu          *\n";
	cout << "* -> mod.dat         zawierajacy najnizszy mod empiryczny sygnalu              *\n";
	cout << "********************************************************************************\n";

	return 0;
}
