// inksci.com，墨之科技

#define G 9.81 // 引力系数
#define ToRadian  0.0174533 // 角度转弧度
#define ToDegree 57.2957805 // 弧度转角度
#define RUN_PACE 0.1 // simulation time step: 0.1 s.
#define NumbVars 10

void main_missile();
void read_data(void);
void rk(double);

extern double y[NumbVars],y00[NumbVars];

