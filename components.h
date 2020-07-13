float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));

    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));

    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}

void calculateZeta(int i, mesh m, Matrix &Ze){
    zeroes(Ze, 12, 3);

    element e = m.getElement(i);

    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);
    //transponiendo x,y
    float x1, x2, x3, x4, y1, y2, y3, y4;
    x1= n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    y1= n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    //Definiendo cada elemento de matriz Z

    //primera columna
    Ze.at(0).at(0) = ((2*x1)+x2+x3+x4-(2*y1)-y2-y3-y4);
    Ze.at(1).at(0) = (x1+(2*x2)+x3+x4-y1-(2*y2)-y3-y4);
    Ze.at(2).at(0) = (x1+x2+(2*x3)+x4-y1-y2-(2*y3)-y4);
    Ze.at(3).at(0) = (x1+x2+x3+(2*x4)-y1-y2-y3-(2*y4));
    Ze.at(4).at(0) = 0;
    Ze.at(5).at(0) = 0;
    Ze.at(6).at(0) = 0;
    Ze.at(7).at(0) = 0;
    Ze.at(8).at(0) = 0;
    Ze.at(9).at(0) = 0;
    Ze.at(10).at(0) = 0;
    Ze.at(11).at(0) = 0;

    //Segunda columna
    Ze.at(0).at(1) =0;
    Ze.at(1).at(1) =0;
    Ze.at(2).at(1) =0;
    Ze.at(3).at(1) =0;
    Ze.at(4).at(0) = ((2*x1)+x2+x3+x4-(2*y1)-y2-y3-y4);
    Ze.at(5).at(0) = (x1+(2*x2)+x3+x4-y1-(2*y2)-y3-y4);
    Ze.at(6).at(0) = (x1+x2+(2*x3)+x4-y1-y2-(2*y3)-y4);
    Ze.at(7).at(0) = (x1+x2+x3+(2*x4)-y1-y2-y3-(2*y4));
    Ze.at(8).at(1) =0;
    Ze.at(9).at(1) =0;
    Ze.at(10).at(1) =0;
    Ze.at(11).at(1) =0;

    //Tercera columna
    Ze.at(0).at(2) =0;
    Ze.at(1).at(2) =0;
    Ze.at(2).at(2) =0;
    Ze.at(3).at(2) =0;
    Ze.at(4).at(2) =0;
    Ze.at(5).at(2) =0;
    Ze.at(6).at(2) =0;
    Ze.at(7).at(2) =0;
    Ze.at(8).at(2) = ((2*x1)+x2+x3+x4-(2*y1)-y2-y3-y4);
    Ze.at(9).at(2) = (x1+(2*x2)+x3+x4-y1-(2*y2)-y3-y4);
    Ze.at(10).at(2) = (x1+x2+(2*x3)+x4-y1-y2-(2*y3)-y4);
    Ze.at(11).at(2) = (x1+x2+x3+(2*x4)-y1-y2-y3-(2*y4));

}

void calculateW(int i,mesh m, Matrix &w){
    zeroes(w, 3, 12);

    element e = m.getElement(i);

    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);
    //transponiendo x
    float x1, x2, x3, x4;
    x1= n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    //primera fila
    w.at(0).at(0) = -(x1 + x2 + x3 + x4);
    w.at(0).at(1) = x1 + x2 + x3 + x4;
    w.at(0).at(2) = 0;
    w.at(0).at(3) = 0;
    w.at(0).at(4) = -(x1 + x2 + x3 + x4);
    w.at(0).at(5) = x1 + x2 + x3 + x4;
    w.at(0).at(6) = 0;
    w.at(0).at(7) = 0;
    w.at(0).at(8) = -(x1 + x2 + x3 + x4);
    w.at(0).at(9) = x1 + x2 + x3 + x4;
    w.at(0).at(10) = 0;
    w.at(0).at(11) = 0;

    //Segunda fila
    w.at(1).at(0) =  -(x1 + x2 + x3 + x4);
    w.at(1).at(1) =  0;
    w.at(1).at(2) =  (x1 + x2 + x3 + x4);
    w.at(1).at(3) =  0;
    w.at(1).at(4) =  -(x1 + x2 + x3 + x4);
    w.at(1).at(5) =  0;
    w.at(1).at(6) =  (x1 + x2 + x3 + x4);
    w.at(1).at(7) =  0;
    w.at(1).at(8) =  -(x1 + x2 + x3 + x4);
    w.at(1).at(9) =  0;
    w.at(1).at(10) =  (x1 + x2 + x3 + x4);
    w.at(1).at(11) =  0;

    //Tercera fila
    w.at(2).at(0) = -(x1 + x2 + x3 + x4);
    w.at(2).at(1) = 0;
    w.at(2).at(2) = 0;
    w.at(2).at(3) = (x1 + x2 + x3 + x4);
    w.at(2).at(4) = -(x1 + x2 + x3 + x4);
    w.at(2).at(5) = 0;
    w.at(2).at(6) = 0;
    w.at(2).at(7) = (x1 + x2 + x3 + x4);
    w.at(2).at(8) = -(x1 + x2 + x3 + x4);
    w.at(2).at(9) = 0;
    w.at(2).at(10) = 0;
    w.at(2).at(11) = (x1 + x2 + x3 + x4);

}

void calculateL(int i, mesh m,Matrix &L){
    zeroes(L, 12, 3);

    element e = m.getElement(i);

    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    //Transponiendo x
    float y1,y2,y3,y4;
    y1= n1.getY();
    y2 = n2.getY();
    y3 = n3.getY();
    y4 = n4.getY();

    //Primera columna

    L.at(0).at(0) = ((2*y1)+y2+y3+y4);
    L.at(1).at(0) = (y1+(2*y2)+y3+y4);
    L.at(2).at(0) = (y1+y2+(2*y3)+y4);
    L.at(3).at(0) = (y1+y2+y3+(2*y4));
    L.at(4).at(0) = 0;
    L.at(5).at(0) = 0;
    L.at(6).at(0) = 0;
    L.at(7).at(0) = 0;
    L.at(8).at(0) = 0;
    L.at(9).at(0) = 0;
    L.at(10).at(0) = 0;
    L.at(11).at(0) = 0;

    //segunda columna

    L.at(0).at(1) = 0;
    L.at(1).at(1) = 0;
    L.at(2).at(1) = 0;
    L.at(3).at(1) = 0;
    L.at(4).at(1) = ((2*y1)+y2+y3+y4);
    L.at(5).at(1) = (y1+(2*y2)+y3+y4);
    L.at(6).at(1) = (y1+y2+(2*y3)+y4);
    L.at(7).at(1) = (y1+y2+y3+(2*y4));
    L.at(8).at(1) = 0;
    L.at(9).at(1) = 0;
    L.at(10).at(1) = 0;
    L.at(11).at(1) = 0;

    //tercera columna
    L.at(0).at(1) = 0;
    L.at(1).at(1) = 0;
    L.at(2).at(1) = 0;
    L.at(3).at(1) = 0;
    L.at(4).at(0) = 0;
    L.at(5).at(0) = 0;
    L.at(6).at(0) = 0;
    L.at(7).at(0) = 0;
    L.at(8).at(1) = ((2*y1)+y2+y3+y4);
    L.at(9).at(1) = (y1+(2*y2)+y3+y4);
    L.at(10).at(1) = (y1+y2+(2*y3)+y4);
    L.at(11).at(1) = (y1+y2+y3+(2*y4));


}

void calculateS(int i, mesh m, Matrix &S){
    zeroes(S, 3, 12);

    element e = m.getElement(i);

    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    //Transponiendo x
    float x1,x2,x3,x4;
    x1= n1.getX();
    x2 = n2.getX();
    x3 = n3.getX();
    x4 = n4.getX();

    //primera fila
    S.at(0).at(0) = ((3*pow(x1,2))+(2*x1*(x2+x3+x4))+pow(x2,2)+(x2*(x3+x4))+pow(x3,2)+(x3*x4)+pow(x4,2));
    S.at(0).at(1) = (pow(x1,2)+(x1*(2*x2+x3+x4))+(3*pow(x2,2))+(2*x2*(x3+x4))+pow(x3,2)+x3*x4+pow(x4,2));
    S.at(0).at(2) = (pow(x1,2)+(x1*(x2+x3+x4))+pow(x2,2)+(x2*(2*x3+x4))+(3*pow(x3,2))+(2*x3*x4)+pow(x4,2));
    S.at(0).at(3) = (pow(x1,2)+(x1*(x2+x3+x4))+pow(x2,2)+(x2*(x3+x4))+pow(x3,2)+(2*x3*x4)+(3*pow(x4,2)));
    S.at(0).at(4) = 0;
    S.at(0).at(5) = 0;
    S.at(0).at(6) = 0;
    S.at(0).at(7) = 0;
    S.at(0).at(8) = 0;
    S.at(0).at(9) = 0;
    S.at(0).at(10) = 0;
    S.at(0).at(11) = 0;

    //segunda fila
    S.at(1).at(0) = 0;
    S.at(1).at(1) = 0;
    S.at(1).at(2) = 0;
    S.at(1).at(3) = 0;
    S.at(1).at(4) = ((3*pow(x1,2))+(2*x1*(x2+x3+x4))+pow(x2,2)+(x2*(x3+x4))+pow(x3,2)+(x3*x4)+pow(x4,2));
    S.at(1).at(5) = (pow(x1,2)+(x1*(2*x2+x3+x4))+(3*pow(x2,2))+(2*x2*(x3+x4))+pow(x3,2)+x3*x4+pow(x4,2));
    S.at(1).at(6) = (pow(x1,2)+(x1*(x2+x3+x4))+pow(x2,2)+(x2*(2*x3+x4))+(3*pow(x3,2))+(2*x3*x4)+pow(x4,2));
    S.at(1).at(7) = (pow(x1,2)+(x1*(x2+x3+x4))+pow(x2,2)+(x2*(x3+x4))+pow(x3,2)+(2*x3*x4)+(3*pow(x4,2)));
    S.at(1).at(8) = 0;
    S.at(1).at(9) = 0;
    S.at(1).at(10) = 0;
    S.at(1).at(11) = 0;
    //tercera fila
    S.at(2).at(0) = 0;
    S.at(2).at(1) = 0;
    S.at(2).at(2) = 0;
    S.at(2).at(3) = 0;
    S.at(2).at(4) = 0;
    S.at(2).at(5) = 0;
    S.at(2).at(6) = 0;
    S.at(2).at(7) = 0;
    S.at(2).at(8) = ((3*pow(x1,2))+(2*x1*(x2+x3+x4))+pow(x2,2)+(x2*(x3+x4))+pow(x3,2)+(x3*x4)+pow(x4,2));
    S.at(2).at(9) = (pow(x1,2)+(x1*(2*x2+x3+x4))+(3*pow(x2,2))+(2*x2*(x3+x4))+pow(x3,2)+x3*x4+pow(x4,2));
    S.at(2).at(10) = (pow(x1,2)+(x1*(x2+x3+x4))+pow(x2,2)+(x2*(2*x3+x4))+(3*pow(x3,2))+(2*x3*x4)+pow(x4,2));
    S.at(2).at(11) = (pow(x1,2)+(x1*(x2+x3+x4))+pow(x2,2)+(x2*(x3+x4))+pow(x3,2)+(2*x3*x4)+(3*pow(x4,2)));

}

float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){

    Matrix matrixQ,matrixH,matrixR,matrixT;
    float u_bar,nu,rho,Ve,J,Determinant;
    
    /* [ Q+H  R ]
       [  T   0 ]
    */

    //Matrix Q
    Matrix ze, Alpha, Beta;
    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);

    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(120*Determinant);

    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    calculateZeta(e,m,ze);
    productRealMatrix(real_a, productMatrixMatrix(ze,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixQ);


    //Matrix H
    Matrix Alpha_t,w_t, doblew, Beta_t;
    
    float real_k = (float) (J)/(24 *Determinant*Determinant);
    calculateW(e,m,doblew);
    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);
    transpose(doblew,w_t);

    productRealMatrix(real_k,productMatrixMatrix(w_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixH);

    
    //Matrix R
    Matrix ele, Omega;
    calculateOmega(Omega);

    float real_g = (float) (J/(120*Determinant));

    calculateL(e,m,ele);
    productRealMatrix(real_g,productMatrixMatrix(ele,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixR);

    //Matrix T
    Matrix S,Omega_t, s_t;
    float real_d = (float)(J/(360*Determinant));
    calculateS(e,m, S);
    transpose(Omega, Omega_t);
    transpose(S, s_t);
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,S,3,3,12),4,3,12),matrixT);

    //Matrix M
    Matrix M;
    zeroes(M,16);
    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixQ,matrixH,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixR);
    ubicarSubMatriz(M,12,15,0,11,matrixT);

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}

Vector createLocalb(int e,mesh &m){
    float J;
    Vector b,b_aux,f;
    Matrix g_matrix;

    calculateF(f, m);

    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    zeroes(b_aux,16);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    
    return b;
}
