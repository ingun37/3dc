#define screenwidth 160
#define screenheight 80
#define pow2(d) (d*d)
enum
{
	vposition=1,
	vdirection=0,
};
enum
{
	cix=0,
	ciy=1,
	ciz=2,
	ciw=3,
};
typedef double* vector;
typedef double** matrix;

#include <stdio.h>
#include <Windows.h>
#include <math.h>

static char screen[screenheight][screenwidth];
static bool overed[screenheight][screenwidth];
vector makeVector(double x, double y, double z, double w)
{
	vector ret = (double*)malloc(sizeof(double)*4);
	ret[0]=x;
	ret[1]=y;
	ret[2]=z;
	ret[3]=w;
	return ret;
}


void freevec(vector vec)
{
	free(vec);
}
void copyvec(vector src, vector dst)
{
	int i;
	for(i=0;i<4;i++)
		dst[i]=src[i];
}
void freemat(matrix mat)
{
	int i;
	for(i=0;i<4;i++)
	{
		free(mat[i]);
	}
	free(mat);
}
matrix makeMatrix(double m11, double m12, double m13, double m14,
				  double m21, double m22, double m23, double m24,
				  double m31, double m32, double m33, double m34,
				  double m41, double m42, double m43, double m44)
{
	int i=0;
	matrix m;
	m = (double**)malloc(sizeof(void*)*4);
	for(i=0;i<4;i++)
	{
		m[i] = (double*)malloc(sizeof(double)*4);
	}
	m[0][0] = m11; m[0][1] = m12; m[0][2] = m13; m[0][3] = m14;
	m[1][0] = m21; m[1][1] = m22; m[1][2] = m23; m[1][3] = m24;
	m[2][0] = m31; m[2][1] = m32; m[2][2] = m33; m[2][3] = m34;
	m[3][0] = m41; m[3][1] = m42; m[3][2] = m43; m[3][3] = m44;
	return m;
}

matrix makeTranslateMatrix(double dx, double dy, double dz)
{
	return makeMatrix(	1,0,0,dx,
						0,1,0,dy,
						0,0,1,dz,
						0,0,0,1);
}

matrix makequaternionrotationmatrix( vector axis, double angle)
{
	double sinangle = sin(angle);
	double cosangle = cos(angle);

	vector q = makeVector( axis[0] * sinangle, axis[1] * sinangle, axis[2] * sinangle, cosangle);
	
	matrix ret = makeMatrix(	1-2*(pow2(q[ciy]) + pow2(q[ciz])), 2*(q[cix]*q[ciy] - q[ciw]*q[ciz]), 2*(q[cix]*q[ciz] + q[ciw]*q[ciy]), 0,
						2*(q[cix]*q[ciy] + q[ciw]*q[ciz]), 1-2*(pow2(q[cix]) + pow2(q[ciz])), 2*(q[ciy]*q[ciz] - q[ciw]*q[cix]), 0,
						2*(q[cix]*q[ciz] - q[ciw]*q[ciy]), 2*(q[ciy]*q[ciz] + q[ciw]*q[cix]), 1 - 2*(pow2(q[cix]) + pow2(q[ciy])), 0,
						0,0,0,1);

	freevec(q);
	return ret;
}

vector mulMatVec(matrix mat,vector vec)
{
	int i, j;
	vector retvec = makeVector(0,0,0,0);
	
	for(i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			retvec[i] += vec[j]*mat[i][j];
		}
	}
	return retvec;
}

vector rotateVector(vector v, vector axis, double angle)
{
	matrix qmat = makequaternionrotationmatrix(axis, angle);
	vector res = mulMatVec(qmat,v);
	freemat(qmat);
	return res;
}

matrix mulMatMat(matrix mat1, matrix mat2)
{
	matrix retmat = makeMatrix(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	int i,j,k;
	for(i=0;i<4;i++)
	{
		for(j=0;j<4;j++)
		{
			for(k=0;k<4;k++)
			{
				retmat[i][j]+=mat1[i][k]*mat2[k][j];
			}
		}
	}
	return retmat;
}

vector cross(vector v1, vector v2)
{
	return makeVector(	v1[ciy]*v2[ciz]-v1[ciz]*v2[ciy],
						v1[ciz]*v2[cix]-v1[cix]*v2[ciz],
						v1[cix]*v2[ciy]-v1[ciy]*v2[cix],
						vdirection);
}

double inner(vector v1, vector v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void makescreen()
{
	int i, j;
	
	HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	COORD coord = {0, 0};
	DWORD count;
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	GetConsoleScreenBufferInfo(hStdOut, &csbi);
	FillConsoleOutputCharacter(hStdOut, ' ',
                               csbi.dwSize.X * csbi.dwSize.Y,
                               coord, &count);
	SetConsoleCursorPosition(hStdOut, coord);

	for(i=0;i<screenheight;i++)
	{
		for(j=0;j<screenwidth;j++)
		{
			printf("%c",screen[i][j]);
		}
		putchar('\n');
	}
}

bool isInFace(vector* startp, unsigned int x, unsigned int y)
{
	vector v1, v2;
	vector v2v, v2p;
	int i;
	double crossres;

	int max=1,min=1;
	for(i=0;i<3;i++)
	{
		if(startp[i][0] >= x)
			max = 0;
		else
			min = 0;
	}
	if(max || min)
		return false;
	for(i=0;i<3;i++)
	{
		v1 = startp[i];
		v2 = startp[(i+1)%3];

		v2v = makeVector(v2[0] - v1[0], v2[1] - v1[1],0,vdirection);
		v2p = makeVector(x-v1[0], y-v1[1], 0, vdirection);
		

		crossres = v2p[0]*v2v[1] - v2v[0]*v2p[1];
		if(crossres > 0)
			return false;//outter point


		freevec(v2v);
		freevec(v2p);
	}

	return true;
}
vector* worldworld(vector* vectors, unsigned int vectornum, matrix translate, matrix rotate, matrix scale)
{
	vector tmpv;
	vector* ret = (vector*)malloc(sizeof(vector)*vectornum);
	int i;
	for(i=0;i<vectornum;i++)
	{
		if(scale)
			ret[i] = mulMatVec(scale,vectors[i]);
		else
			ret[i] = makeVector(vectors[i][0],vectors[i][1],vectors[i][2],vectors[i][3]);

		if(rotate)
		{
			tmpv = ret[i];
			ret[i] = mulMatVec(rotate, tmpv);
			freevec(tmpv);
		}

		if(translate)
		{
			tmpv = ret[i];
			ret[i] = mulMatVec(translate,tmpv);
			freevec(tmpv);
		}

	}
	return ret;
}
vector* worldview(vector* vectors, unsigned int vectornum, vector campos, vector camF)
{
	int i;
	vector tmpv, camR, camU;
	vector aUp;

	vector* ret = (vector*)malloc(sizeof(vector)*vectornum);
	aUp = makeVector(0,1,0,vdirection);
	camR = cross(camF,aUp);
	camU = cross(camR,camF);
	freevec(aUp);

	matrix rotmat = makeMatrix(	camR[0],camR[1],camR[2],0,
								camU[0],camU[1],camU[2],0,
								camF[0],camF[1],camF[2],0,
								0,0,0,1);
	matrix movmat = makeMatrix(	1,0,0,-campos[0],
								0,1,0,-campos[1],
								0,0,1,-campos[2],
								0,0,0,1);

	matrix viewmat = mulMatMat(rotmat,movmat);
	freemat(rotmat);
	freemat(movmat);

	for(i=0;i<vectornum;i++)
	{
		tmpv = vectors[i];
		ret[i] = mulMatVec(viewmat,vectors[i]);
		//vectors[i][ciz]=-vectors[i][ciz];
		//freevec(tmpv);
	}

	return ret;
}

vector* worldproj(vector* vectors, unsigned int vectornum, double r, double l, double b, double t, double n, double f)
{
	int i;
	vector tmpv;
	vector* ret = (vector*)malloc(sizeof(vector)*vectornum);
	matrix proj = makeMatrix(	3.f/2,0,0,0,
								0,1,0,0,
								0,0,1,0,
								0,0,-1/n,0);
	for(i=0;i<vectornum;i++)
	{
		tmpv = vectors[i];
		ret[i] = mulMatVec(proj,tmpv);
		//freevec(tmpv);
	}

	return ret;
}

vector* worldscr(vector* vectors, unsigned int vectornum, double r, double l, double t, double b)
{
	int i,j;
	double w;
	vector tmpv;
	vector* ret = (vector*)malloc(sizeof(vector)*vectornum);
	matrix scr = makeMatrix(	(screenwidth/2)/r,0,0,(screenwidth/2),
								0,-(screenheight/2)/t,0,(screenheight/2),
								0,0,1,0,
								0,0,0,1);

	for(i=0;i<vectornum;i++)
	{
		tmpv = makeVector( vectors[i][0],vectors[i][1],vectors[i][2],vectors[i][3] );
		w = tmpv[ciw];
		for(j=0;j<4;j++)
		{
			tmpv[j] = tmpv[j]/w;
		}

		ret[i] = mulMatVec(scr,tmpv);
		freevec(tmpv);
	}

	return ret;
}
int main()
{
	int i, j, rs, k;
	unsigned int facenum = 1;
	char input;
	matrix trans, rot;
	double angle = 0;
	double angle2 = 0;
	const double startingxpos = -15;
	double xpos = startingxpos;
	vector axis;
	vector* vectors;
	vector *afterworld, *afterview, *afterproj, *afterscr, *startp;
	vector campos;
	vector camf;
	int stage=0;
	int stageflag[2]={0,0};
	double tmpd;
	for(i=0;i<screenheight;i++)
		for(j=0;j<screenwidth;j++)
		{
			screen[i][j]='_';
			overed[i][j] = false;
		}

	campos = makeVector(0,0,-1,vposition);
	camf = makeVector(0,0,-1,vdirection);

	
	
	while(1)
	{
		switch(stage)
		{
		case 0:
			if(stageflag[stage]==0)
			{
				facenum = 1;

				vectors = (vector*)malloc(sizeof(vector)*3*facenum);
	
				tmpd = sqrt(3.f);
				vectors[0] = makeVector(-2,-2*tmpd*(1.f/3),0,vposition);
				vectors[1] = makeVector(0,2*tmpd*(2.f/3),0,vposition);
				vectors[2] = makeVector(2,-2*tmpd*(1.f/3),0,vposition);

				stageflag[stage]=1;
			}
			trans = makeTranslateMatrix(xpos,0,-10);
			axis = makeVector(0,0,1,vdirection);
			rot = makequaternionrotationmatrix(axis,angle);
			freevec(axis);

			afterworld = worldworld(vectors,facenum*3,trans,rot,NULL);
			break;
		case 1:
			if(stageflag[stage]==0)
			{
				for(i=0;i<facenum*3;i++)
				{
					freevec(vectors[i]);
				}

				free(vectors);

				facenum=3;

				vectors = (vector*)malloc(sizeof(vector)*3*facenum);
	
				tmpd = sqrt(3.f);
				vectors[0] = makeVector(-2,-2*tmpd*(1.f/3) + 2*tmpd*(2.f/3),0,vposition);
				vectors[1] = makeVector(0,2*tmpd*(2.f/3) + 2*tmpd*(2.f/3),0,vposition);
				vectors[2] = makeVector(2,-2*tmpd*(1.f/3) + 2*tmpd*(2.f/3),0,vposition);

				vectors[3] = makeVector(-2 - 2,-2*tmpd*(1.f/3) - 2*tmpd*(1.f/3),0,vposition);
				vectors[4] = makeVector(0 - 2,2*tmpd*(2.f/3) - 2*tmpd*(1.f/3),0,vposition);
				vectors[5] = makeVector(2 - 2,-2*tmpd*(1.f/3) - 2*tmpd*(1.f/3),0,vposition);

				vectors[6] = makeVector(-2 + 2,-2*tmpd*(1.f/3) - 2*tmpd*(1.f/3),0,vposition);
				vectors[7] = makeVector(0 + 2,2*tmpd*(2.f/3) - 2*tmpd*(1.f/3),0,vposition);
				vectors[8] = makeVector(2 + 2,-2*tmpd*(1.f/3) - 2*tmpd*(1.f/3),0,vposition);

				trans = makeTranslateMatrix(0,0,-50);

				stageflag[stage]=1;
			}
			axis = makeVector(0,1,0,vdirection);
			rot = makequaternionrotationmatrix(axis,angle2);

			afterworld = worldworld(vectors,facenum*3,trans,rot,NULL);
			break;
		}
		
		afterview = worldview(afterworld,facenum*3,campos,camf);
		afterproj = worldproj(afterview,facenum*3,2,-2,-1,1,-1,-1000);
		afterscr = worldscr(afterproj,facenum*3,2,-2,1,-1);

		for(i=0;i<screenheight;i++)
		{
				for(j=0;j<screenwidth;j++)
				{
					if(overed[i][j])
						screen[i][j] = ' ';
					else
						screen[i][j] = '_';

					for(k = 0;k<facenum;k++)
					{
						startp = afterscr + (k*3);
						if(isInFace(startp,j,i))
						{
							screen[i][j] = 'a' + rand()%('Z'-'a');
							overed[i][j] = true;
						}
					}
					
				}
		}
		makescreen();
		/*
		if(afterscr[0][0] >= 0 && afterscr[0][1] >= 0 && afterscr[0][0] < screenwidth && afterscr[0][1] < screenheight)
		{
			screen[(int)afterscr[0][1]][(int)afterscr[0][0]] = 'X';
			makescreen();
		}
		*/
		input = getchar();
		switch(input)
		{
		case 'd':
			campos[0]+=0.3f;
			break;
		case 'a':
			campos[0]-=0.3f;
			break;
		case 'w':
			campos[1]+=0.3f;
		break;
		case 's':
			campos[1]-=0.3f;
		break;
		default:
			switch(stage)
			{
			case 0:
				angle+=0.1;
				xpos += 0.5f;
				if(xpos >= -startingxpos)
					stage=1;
				break;
			case 1:
				campos[ciz]-=0.8f;
				angle2+=0.05;
				break;
			}
			break;
		}
	}
	return 0;
}