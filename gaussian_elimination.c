// Написать программу, которая будет считывать из файла коэффициенты матрицы и решать её методом Гаусса 

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#define min 0.0001

void fillMatrix(FILE* f, double** m, int sz1, int sz2); 
// Функция, считывающая значения коэффициентов из файла и заполяющая ими матрицу
double* create1(int sz); // Функция выделения памяти под массив решения системы
double** create2(int sz1, int sz2); // Функция, создающая двумерный массив
void erase2(double** m, int sz1, int sz2); // Функция, удаляющая из памяти двумерный массив
void print(double** m, int sz1); // Функция печати двумерного массива
int triangleMatrix(double** m, int sz1); // Функция, приводящая СУ к треугольному виду
void computeAnswer(double** m, double* x, int sz); // Функция, реализующая обратную подстановку
int findLine(double** m, int sz1, int k); // Функция, проверяющая невырожденность системы 
void swap(double** m, int b, int k); // Функция, меняющая строки местами
void countKoeff(int k, int i, int sz1, double n, double** m); // Пересчёт коэффициентов
void printResult(double** m, double* x, int sz1, int val); // Вывод результата вычислений
int getMatrixType(double** m, int sz1); // Функция, определяющая причину, по которой матрица не имеет решений

int getMatrixType(double** m, int sz1) // Функция, определяющая причину, по которой матрица не имеет решений
{
	int i, j;
	for(i=1; i<sz1; i++) // Начиная со 2-го уравнения системы, проверяем каждый коэффициент:
		for (j=0; j<sz1+1; j++)
		{
			if (m[i][j]!=0)  // Если встретилось число, отличное от нуля, проверяем, где оно стоит
				if (j==sz1)
					return 3; // Если это свободный член, значит, у системы нет ни одного решения
				else 
					return 0; // Если это коэффициент при переменной, значит, система вырождена
		} 
	return 2; // Если все нули, значит, коэффициенты системы пропорциональны, и она имеет бесконечное кол-во решений
}

void printResult(double** m, double* x, int sz1, int val) // Вывод результата вычислений
{
	switch (val){
		case 0: 
		{
			printf("Singular matrix"); // Система вырождена
			break;
		}
		case 1: 
		{
			computeAnswer(m, x, sz1); // Обратная подстановка
			int i;
			for(i=0; i<sz1; i++)
				printf("x%i=%15E\n", i, x[i]);
			break;
		}
		case 2: 
		{
			printf("System has endless amount of solutions"); // Система имеет бесконечное множество решений
			break;
		}
		case 3: 
		{
			printf("System has no solutions"); // Система не имеет ни одного решения
			break;
		}
	}
}

void countKoeff(int k, int i, int sz1, double n, double** m) // Пересчёт коэффициентов
{
	int j;
	double x;
	for (j=k+1; j<sz1+1; j++){
		x=m[i][j]-n*m[k][j]; 
		m[i][j]=(fabs(x)<min)?0:x;
	}
}

void swap(double** m, int b, int k) // Функция, меняющая строки местами
{
	double* q=m[k];
	m[k]=m[b];
	m[b]=q;	
}

int findLine(double** m, int sz1, int k) // Функция, которая ищет строку с первым ненулевым коэфф. для обмена
{
	int i;
	for (i=k; i<sz1; i++) // Для каждого уравнения, начиная с текущего, выполняется проверка
	{
		if (fabs(m[i][k])>min)
			return i; // Если хотя бы один коэффициент, стоящий первым, не равен 0, значит, возвращаем № этой строки
	}
	return 0; // Возвращается 0, если система вырождена
}

int triangleMatrix(double** m, int sz1) // Функция, приводящая СУ к треугольному виду
{ 
	int i, k, j, b;
	double n;
	for (k=0; k<sz1; k++) // Выполняем действия для каждого уравнения:
	{ 
		for (i=k+1; i<sz1; i++) // Начиная со 2-го уравнения системы, исключаем неизвестные из уравнений
		{
			if (fabs(m[k][k])<min) // Если в текущей строке первый коэффициент равен 0, то ищем строку для обмена
			{
				if (b=findLine(m, sz1, k))
					swap(m, b, k); // Если строка с 1ым ненулевым коэффициентом нашлась, то меняем строки местами
			 	else 
			 		return getMatrixType(m, sz1);
			}
			n=m[i][k]/m[k][k]; // Вводится коэффициент
			m[i][k]=0; // 1-ый элемент в строке исключён, поэтому присваиваем 0
			countKoeff(k, i, sz1, n, m);
		}
		print(m, sz1);
	}	
	return 1; // Система имеет решение
}

void computeAnswer(double** m, double* x, int sz1) // Функция, реализующая обратную подстановку
{ 
	int k, j;
	double sum, n;
	for (j=sz1-1; j>=0; j--) // Выражаем корни из каждого уравнения, начиная с конца 
	{
		n=0; 
		for (k=j+1; k<sz1; k++) // Кол-во повторений = кол-во вычитаемых чисел 
		{
			sum=m[j][k]*x[k]; // Корни из следующих уравнений, умноженные на коэффициенты, стоящие перед ними
		//	printf("%15E\n", sum);
			n+=sum; // Собираем все вычитаемые числа
		}
		x[j]=(m[j][sz1]-n)/m[j][j]; // Выражаем очередной корень уравнения
	}
}

void fillMatrix(FILE* f, double** m, int sz1, int sz2)
{ // Функция, считывающая значения коэффициентов из файла и заполняющая ими матрицу
	int i, j;
	for (i=0; i<sz1; i++)
		for (j=0; j<sz2; j++)
			fscanf(f, "%lf", &m[i][j]);
}

double* create1(int sz) // Функция выделения памяти под массив решения системы
{
	double* m=(double*) malloc (sz*sizeof (double));
	return m;
}

double** create2(int sz1, int sz2) // Функция, создающая двумерный массив
{ 
	double** m=(double**) malloc (sz1*sizeof (double*)); // Выделяется память под массив указателей
	int i;
	for (i=0; i<sz1; i++) 
		m[i]=(double*) malloc (sz2*sizeof (double)); // Выделяется память под массив элементов типа double
	return m;
}

void erase2(double** m, int sz1, int sz2) // Функция, удаляющая из памяти двумерный массив
{ 
	int i;
	for (i=0; i<sz1; i++) 
		free(m[i]); // Очищаются массивы элементов
	free(m); // Очищается массив указателей
}

void print(double** m, int sz1) // Функция печати двумерного массива
{ 
	int i, j;
	for (i=0; i<sz1; i++){
		for (j=0; j<sz1+1; j++) 
			printf ("%15E", m[i][j]); 
		printf ("\n");
	}
	printf("\n");
}

int main()
{
	int sz1, sz2, b;
	FILE* f=fopen("gauss-z1.txt", "r"); // Открываем файл для чтения
	assert(f!=NULL); // Проверяем, открылся ли файл
	fscanf(f, "%i %i", &sz1, &sz2); // Считываем размеры матрицы, хранящиеся в 1ой строке файла
	double** m=create2(sz1, sz2); // Выделяем память под матрицу
	double* x=create1(sz1);
	fillMatrix(f, m, sz1, sz2); // Заполняем матрицу значениями из файла
	fclose(f);
	print(m, sz1);
	b=triangleMatrix(m, sz1); // Если матрица приведена к треугольному виду успешно, то проводятся дальнейшие действия
	print(m, sz1); 
	printResult(m, x, sz1, b);
	erase2(m, sz1, sz2);
	free(x);
	return 0;
}
