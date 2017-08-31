/// 雅可比方法求实对称矩阵的特征值和特征向量
/// 
/// \param[in]   n       number of dimension
/// \param[in]   eps     precision
/// \param[in]   jt      maximum iteration number
/// \param[in,out]   a   initial matrix / diagonal eigenvalue matrix
/// \param[out]  v       column eigenvector
int cjcbi(real a[][], int n, real v[][], real eps, int jt) {
    int i,j,p,q,u,w,t,s,l;
    real fm,cn,sn,omega,x,y,d;
    l=1;
    for(i = 0; i <= n-1; ++i) {
        v[i][i] = 1.0;
        for(j = 0; j <= n-1; ++j)
          if (i != j) v[i][j] = 0.0;
    }
    while (1 == 1) {
        fm = 0.0;
        for(i = 1; i <= n-1; ++i) {
            for(j = 0; j <= i-1; ++j) {
                d = fabs(a[i][j]);
                if((i != j)&&(d > fm)) {
                    fm = d; p = i; q = j;
                }
            }
        }
        if(fm < eps)  return(1);
        if(l > jt)  return(-1);
        l = l+1;
        u = p*n+q; w = p*n+p; t = q*n+p; s = q*n+q;
        x = -a[p][q]; y = (a[q][q]-a[p][p])/2.0;
        omega = x/sqrt(x*x+y*y);
        if(y < 0.0) omega = -omega;
        sn = 1.0 + sqrt(1.0 - omega*omega);
        sn = omega/sqrt(2.0*sn);
        cn = sqrt(1.0 - sn*sn);
        fm = a[p][p];
        a[p][p] = fm*cn*cn + a[q][q]*sn*sn + a[p][q]*omega;
        a[q][q] = fm*sn*sn + a[q][q]*cn*cn - a[p][q]*omega;
        a[p][q] = 0.0; a[q][p] = 0.0;
        for(j = 0; j <= n-1; ++j) {
            if((j != p) && (j != q)) {
                u =p*n+j; w=q*n+j;
                fm = a[p][j];
                a[p][j] =  fm*cn + a[q][j]*sn;
                a[q][j] = -fm*sn + a[q][j]*cn;
            }
        }
        for(i = 0; i <= n-1; ++i) {
            if((i != p) && (i != q)) {
                u=i*n+p; w=i*n+q;
                fm=a[i][p];
                a[i][p] =  fm*cn + a[i][q]*sn;
                a[i][q] = -fm*sn + a[i][q]*cn;
            }
        }
        for(i = 0; i <= n-1; ++i) {
            u=i*n+p; w=i*n+q;
            fm=v[i][p];
            v[i][p] =  fm*cn + v[i][q]*sn;
            v[i][q] = -fm*sn + v[i][q]*cn;
        }
    }
    return(1);
}

/// Produce the ellipse path based on mean value and covariance matrix
/// 
/// \param[in]  mean        mean value vector
/// \param[in]  Cov         covariance matrix
/// \param[in]  Xindex      index on x-axis
/// \param[in]  Yindex      index on Y-axis
/// \return     the path of ellipse
path Solve_Ellipse(real mean[], real Cov[][], int Xindex, int Yindex) {
    path ellip;

    int jt = 100;
    real eps = 1e-6;

    real Mtx_2D[][] = new real [2][2],
         Vec_2D[][] = new real [2][2];

    pair cent = (mean[Xindex], mean[Yindex]);

    Mtx_2D[0][0] = Cov[Xindex][Xindex];
    Mtx_2D[0][1] = Cov[Xindex][Yindex];
    Mtx_2D[1][0] = Cov[Yindex][Xindex];
    Mtx_2D[1][1] = Cov[Yindex][Yindex];

    // Eigen value and vector
    int flag = cjcbi(Mtx_2D, 2, Vec_2D, eps, jt);
    // ellipse
    if(flag != 0) {
        ellip = rotate(degrees(atan(Vec_2D[1][0]/Vec_2D[0][0])), cent) * 
            ellipse(cent, 3*sqrt(Mtx_2D[0][0]), 3*sqrt(Mtx_2D[1][1]));
    }

    return ellip;
}