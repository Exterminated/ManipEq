#define _USE_MATH_DEFINES
#include "cmath"
#include <conio.h>
#include <iostream>

using namespace std;

void main() {
	double xb, y0, za, xc, yb, tc, zd, O1A, OA1, l1, l2, l3, OO1;
	double  fi, l4, A, gamma1, gamma0, gamma, xm0, ym0, zm0, xe0, ye0, ze0, psi1, psi2, alpha;
	double xmk, ymk, zmk, xek, yek, zek;
	int alpha_0, alpha_23, a, alpha_13, alpha_33, b;
	//�������� �������������� ������
	xb = 355.0;
	y0 = 755.0;
	za = 750.0;
	xc = -355.0;
	yb = y0;
	zd = -40.0;
	O1A = za;
	OA1 = 750.0;
	l1 = 1400.0;
	l2 = 1500.0;
	l3 = 1352.0;
	fi = 0.323;
	OO1 = zd;

	//��������� ���������� ����� � � ����� 4 �����
	l4 = sqrt((za*sin(M_PI_2 - fi) - OO1)*(za*sin(M_PI_2 - fi) - OO1) + (yb - za*cos(M_PI_2 - fi))*(yb - za*cos(M_PI_2 - fi)));
	A = 0.5*l2*l2 + 0.5*l3*l3 - l1*l1;
	xm0 = (l3*l3 - l2*l2) / (4 * xb);	
	ym0 = sqrt((l1*l1) - (((l3*l3 - l2*l2)*(l3*l3 - l2*l2)) / (16 * xb * xb)) - (((A - xb*xb - za*za)*(A - xb*xb - za*za)) / (4 * za*za)))*cos(fi) - ((A - xb*xb + za*za) / (2 * za))*sin(fi);
	zm0 = sqrt((l1*l1) - (((l3*l3 - l2*l2)*(l3*l3 - l2*l2)) / (16 * xb * xb)) - (((A - xb*xb - za*za)*(A - xb*xb - za*za)) / (4 * za*za)))*sin(fi) + ((A - xb*xb + za*za) / (2 * za))*cos(fi);
	gamma0 = atan((-1 * xm0) / (ym0 + OA1*sin(fi)));
	cout << "l4: " << l4 << endl;
	cout << "A: " << A << endl;
	cout << "xm0: " << xm0 << endl;
	cout << "ym0: " << ym0 << endl;
	cout << "zm0: " << zm0 << endl;
	cout << "gamma0: " << gamma0 << endl;
	//��������� ���������� ������� E, ������������ ���������, ��������� �� ���.��������
	alpha_0 = 0;
	alpha_23 = 1;
	a = 90;
	alpha_13 = 0;
	alpha_33 = 0;
	b = 220;

	xe0 = xm0 + b*alpha_13 - a*cos(alpha_0)*sin(gamma0);
	ye0 = ym0 + b*alpha_23 + a*cos(alpha_0)*cos(gamma0);
	ze0 = zm0 + b*alpha_33 + a*sin(alpha_0);

	cout << "xe0: " << xe0 << endl;
	cout << "ye0: " << ye0 << endl;
	cout << "ze0: " << ze0 << endl;
	
	//�������� ���������� ������� �
	gamma = gamma0;
	xek = 300;
	yek = 1500;
	zek = -100;
	
	//���������� ���� psi1, gamma
	psi1 = asin(-alpha_23*sin(gamma) + alpha_13*cos(gamma));
	psi2 = asin(-alpha_23*sin(gamma) - alpha_13*cos(gamma));
	alpha = asin(alpha_33 / cos(psi1));
	gamma1 = atan((-xek) / yek);

	cout << "psi1: " << psi1 << endl;
	cout << "psi2: " << psi2 << endl;
	cout << "alpha: " << alpha << endl;
	cout << "gamma1: " << gamma1 << endl;

	xmk = xek - b*alpha_13 + a*cos(alpha)*sin(gamma1);
	ymk = yek - b*alpha_23 - a*cos(alpha)*cos(gamma1);
	zmk = zek - b*alpha_33 - a*sin(alpha);

	cout << "xmk: " << xmk << endl;
	cout << "ymk: " << ymk << endl;
	cout << "zmk: " << zmk << endl;
	
	//����������� ���� ��
	double OK, OA, OB, DK, Imin;
	OK = y0;
	OA = za;
	OB = xb;
	DK = zd;
	OA1 = za;
	
	Imin = pow((sqrt((xmk*xmk + (ymk + OA1*sin(fi))*(ymk + OA1*sin(fi)) + ((zmk - OA1*cos(fi))*(zmk - OA1*cos(fi))))) - l1), 2.0) + pow((sqrt((OK - OA1*sin(fi))*(OK - OA1*sin(fi)) + (OA1*cos(fi) - DK)*(OA1*cos(fi) - DK)) - l4), 2.0);

	cout << "Imin(fi): " << Imin << endl;

	_getch();
}