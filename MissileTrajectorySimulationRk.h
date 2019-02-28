void rk(double rk_step) {
 double a[4],b[NumbVars],c[NumbVars];
 int i,j;
 a[0]=rk_step/2;a[1]=rk_step/2;a[2]=rk_step;a[3]=rk_step;
 dery(y);
 for(i=0;i<NumbVars;i++) b[i]=y[i];
 for(j=0;j<3;j++)  {
	for(i=0;i<NumbVars;i++) {
		c[i]=b[i]+a[j]*dy[i];
		y[i]=y[i]+a[j+1]*dy[i]/3; }
	dery(c);  }
 for(i=0;i<NumbVars;i++) y[i]=y[i]+a[0]*dy[i]/3; }
