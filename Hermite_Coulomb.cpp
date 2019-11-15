#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <sys/times.h>
#include <vector>
#include <limits>
#include <cfloat>

//#define TREE_STAT 1

using namespace std;

static const int P = 10; // Order of Far-field approximation
static const int Pflat = (P + 1) * (P + 1) * (P + 1);
static const int numpars = 100000; // numpars points total
static const int N0 = 2000;
static const double sq_theta = 0.25; // theta = 0.7
static const double pi = 3.141592653589793238462643383;


double xyzminmax[6];

#ifdef TREE_STAT
int max_leaf_pars = 0;
int min_leaf_pars = 999999999;
double max_leaf_ratio = 1.0;
double min_leaf_ratio = 999999999.0;
double max_leaf_diameter = 0.0;
double min_leaf_diameter = 999999999.0;
int split2_times = 0;
int split4_times = 0;
int split8_times = 0;
int max_level = -1;
int level_cnt[9999] = {0};
#endif

struct xyz // Particle coordinates (physical)
{
  size_t size;
  double* x;
  double* y;
  double* z;
  size_t* index;
  size_t* old_index;

  xyz(size_t numpars_in) : size(numpars_in)
  {
    x = new double[size];
    y = new double[size];
    z = new double[size];
    index = new size_t[size];
    old_index = new size_t[size];
  }

  ~xyz()
  {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] index;
    delete[] old_index;
  }
};

struct panel
{
  size_t members[2];
  double xinterval[2];
  double yinterval[2];
  double zinterval[2];
  double xc; // Panel center x coordinate
  double yc; // Panel center y coordinate
  double zc; // Panel center z coordinate
  vector<size_t> children;
  double MAC; // r^2 / theta^2
  double moments[8][Pflat];
  int moment_flag;
  double t1[P + 1]; // Interpolation points in x direction
  double t2[P + 1];
  double t3[P + 1];
  double wx0[P + 1];
  double wy0[P + 1];
  double wz0[P + 1];

  // Initialization
  panel() : xc(0.0), yc(0.0), zc(0.0), MAC(0.0), moment_flag(0)
  {
    memset(members, 0, sizeof(members));
    memset(xinterval, 0, sizeof(xinterval));
    memset(yinterval, 0, sizeof(yinterval));
    memset(zinterval, 0, sizeof(zinterval));
    memset(moments, 0, sizeof(moments));
    memset(t1, 0, sizeof(t1));
    memset(t2, 0, sizeof(t2));
    memset(t3, 0, sizeof(t3));
    memset(wx0, 0, sizeof(wx0));
    memset(wy0, 0, sizeof(wy0));
    memset(wz0, 0, sizeof(wz0));

    children.reserve(8);
  }
};

struct xyz particles(numpars);

vector<panel> tree;
vector<size_t> leaf;

size_t node_count = 0;

double* lambda;

double* velo;
double* velo_true;

long getTickCount()
{
  tms tm;
  return times(&tm);
}

double minval(double* x)
{
  double MinVal = x[0];
  for (int i = 1; i < numpars; i++) {
    if (MinVal > x[i])
      MinVal = x[i];
  }

  MinVal = MinVal;

  return MinVal;
}

double maxval(double* x)
{
  double MaxVal = x[0];
  for (int i = 1; i < numpars; i++) {
    if (MaxVal < x[i])
      MaxVal = x[i];
  }

  MaxVal = MaxVal;

  return MaxVal;
}

int compute_Ds()
{
  double temp_x, temp_y, temp_z;
  double px, py, pz;
  double xx, yy, zz;
  double sq_xx, sq_yy, sq_zz;
  double sum;
  double R2, R, Rinv;

  for (int i = 0; i < numpars; i++) {
    temp_x = particles.x[i];
    temp_y = particles.y[i];
    temp_z = particles.z[i];

    sum = 0.0;

    for (int j = 0; j < numpars; j++) {
      px = particles.x[j];
      py = particles.y[j];
      pz = particles.z[j];

      xx = temp_x - px;
      yy = temp_y - py;
      zz = temp_z - pz;

      sq_xx = xx * xx;
      sq_yy = yy * yy;
      sq_zz = zz * zz;

      R2 = sq_xx + sq_yy + sq_zz;
      R = sqrt(R2);

      if (i != j) {
        Rinv = 1.0 / R;
        sum += Rinv * lambda[j];
      }
    }

    velo_true[i] = sum;
  }

  return 0;
}

void build_tree_init()
{
  panel temp_panel;

  // Indices of particles belonging to the panel
  temp_panel.members[0] = 0;
  temp_panel.members[1] = numpars - 1;

  // Interval defining the panel
  temp_panel.xinterval[0] = xyzminmax[0];
  temp_panel.xinterval[1] = xyzminmax[1];
  temp_panel.yinterval[0] = xyzminmax[2];
  temp_panel.yinterval[1] = xyzminmax[3];
  temp_panel.zinterval[0] = xyzminmax[4];
  temp_panel.zinterval[1] = xyzminmax[5];

  temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
  temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
  temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);

  double xL = temp_panel.xinterval[1] - temp_panel.xinterval[0];
  double yL = temp_panel.yinterval[1] - temp_panel.yinterval[0];
  double zL = temp_panel.zinterval[1] - temp_panel.zinterval[0];

  double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
  temp_panel.MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

  tree.push_back(temp_panel);
  node_count = 1;
}

void Swap(size_t i, size_t j)
{
  if (i == j)
    return;

  double x = particles.x[i];
  double y = particles.y[i];
  double z = particles.z[i];
  size_t index = particles.index[i];
  size_t old_index = particles.old_index[i];

  double temp_lambda = lambda[i];

  particles.x[i] = particles.x[j];
  particles.y[i] = particles.y[j];
  particles.z[i] = particles.z[j];
  particles.index[i] = particles.index[j];
  particles.old_index[i] = particles.old_index[j];

  lambda[i] = lambda[j];

  particles.x[j] = x;
  particles.y[j] = y;
  particles.z[j] = z;
  particles.index[j] = index;
  particles.old_index[j] = old_index;

  lambda[j] = temp_lambda;
}

void split_2(size_t panel_index, int split_code)
{
  panel child[2];

/*
      -------------------
      |        |        |
      |        |        |
      |Child 0 |Child 1 |
      |        |        |
      |        |        |
-----------------------------> axis A
    start     mid      end
*/

  double tp_x0 = tree[panel_index].xinterval[0];
  double tp_x1 = tree[panel_index].xinterval[1];
  double tp_y0 = tree[panel_index].yinterval[0];
  double tp_y1 = tree[panel_index].yinterval[1];
  double tp_z0 = tree[panel_index].zinterval[0];
  double tp_z1 = tree[panel_index].zinterval[1];

  for (int i = 0; i < 2; i++) {
    child[i].xinterval[0] = tp_x0;
    child[i].xinterval[1] = tp_x1;
    child[i].yinterval[0] = tp_y0;
    child[i].yinterval[1] = tp_y1;
    child[i].zinterval[0] = tp_z0;
    child[i].zinterval[1] = tp_z1;
  }

  double xL = tp_x1 - tp_x0;
  double yL = tp_y1 - tp_y0;
  double zL = tp_z1 - tp_z0;

  double* intervalA[2] = {NULL, NULL};
  double* coordA = NULL;
  double startpointA = 0.0, midpointA = 0.0, endpointA = 0.0;

  if (split_code == 4) { // XYZ = 100, A is X
    xL *= 0.5;
    intervalA[0] = child[0].xinterval;
    intervalA[1] = child[1].xinterval;
    coordA = particles.x;
    startpointA = tp_x0;
    endpointA = tp_x1;
  }
  else if (split_code == 2) { // XYZ = 010, A is Y
    yL *= 0.5;
    intervalA[0] = child[0].yinterval;
    intervalA[1] = child[1].yinterval;
    coordA = particles.y;
    startpointA = tp_y0;
    endpointA = tp_y1;
  }
  else if (split_code == 1) { // XYZ = 001, A is Z
    zL *= 0.5;
    intervalA[0] = child[0].zinterval;
    intervalA[1] = child[1].zinterval;
    coordA = particles.z;
    startpointA = tp_z0;
    endpointA = tp_z1;
  }

  midpointA = 0.5 * (startpointA + endpointA);

  // Child 0 ends with mid point on axis A
  intervalA[0][1] = midpointA;

  // Child 1 begins with mid point on axis A
  intervalA[1][0] = midpointA;

  double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
  double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

  for (int i = 0; i < 2; i++) {
    child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
    child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
    child[i].zc = 0.5 * (child[i].zinterval[0] + child[i].zinterval[1]);
    child[i].MAC = MAC;    
  }
  
  vector<size_t> v[2];
  size_t start = tree[panel_index].members[0];
  size_t end = tree[panel_index].members[1];
  size_t* addr_table = new size_t[end - start + 1];

  size_t index;
  for (index = start; index <= end; index++) {
    particles.index[index] = index;
    addr_table[index - start] = index;

    if (coordA[index] <= midpointA)
      v[0].push_back(index);
    else
      v[1].push_back(index);
  }

  size_t seq = start;
  for (size_t j = 0; j < 2; j++) {
    size_t size = v[j].size();

    if (size >= 1) {
      for (size_t k = 0; k < size; k++) {
        if (k == 0)
          child[j].members[0] = seq;
        if (k == size - 1)
          child[j].members[1] = seq;

        index = v[j][k];
        /*
        // This is very slow
        size_t pos;
        for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++) {
          if (particles.index[pos] == index)
            break;
        }
        Swap(pos, seq);
        */
        // This uses an address table
        size_t pos = addr_table[index - start];
        size_t out = particles.index[seq];
        Swap(pos, seq);
        addr_table[index - start] = seq;
        addr_table[out - start] = pos;

        seq++;
      }

      node_count++;
      tree[panel_index].children.push_back(node_count - 1);
      tree.push_back(child[j]);
      v[j].clear();
    }
  }

  delete[] addr_table;
}

void split_4(size_t panel_index, int XYZ_flag)
{
  panel child[4];

/*
      ^ axis B
      |
  end -------------------
      |        |        |
      |Child 2 |Child 3 |
  mid |------------------
      |        |        |
start |Child 0 |Child 1 |
-----------------------------> axis A
    start     mid      end
*/

  double tp_x0 = tree[panel_index].xinterval[0];
  double tp_x1 = tree[panel_index].xinterval[1];
  double tp_y0 = tree[panel_index].yinterval[0];
  double tp_y1 = tree[panel_index].yinterval[1];
  double tp_z0 = tree[panel_index].zinterval[0];
  double tp_z1 = tree[panel_index].zinterval[1];

  for (int i = 0; i < 4; i++) {
    child[i].xinterval[0] = tp_x0;
    child[i].xinterval[1] = tp_x1;
    child[i].yinterval[0] = tp_y0;
    child[i].yinterval[1] = tp_y1;
    child[i].zinterval[0] = tp_z0;
    child[i].zinterval[1] = tp_z1;
  }

  double xL = tp_x1 - tp_x0;
  double yL = tp_y1 - tp_y0;
  double zL = tp_z1 - tp_z0;

  double* intervalA[4] = {NULL, NULL, NULL, NULL};
  double* intervalB[4] = {NULL, NULL, NULL, NULL};
  double* coordA = NULL;
  double* coordB = NULL;
  double startpointA = 0.0, endpointA = 0.0, midpointA = 0.0;
  double startpointB = 0.0, endpointB = 0.0, midpointB = 0.0;

  switch (XYZ_flag) {
    case 6: // XYZ = 110, A is X and B is Y 
      xL *= 0.5;
      yL *= 0.5;
      intervalA[0] = child[0].xinterval;
      intervalB[0] = child[0].yinterval;
      intervalA[1] = child[1].xinterval;
      intervalB[1] = child[1].yinterval;
      intervalA[2] = child[2].xinterval;
      intervalB[2] = child[2].yinterval;
      intervalA[3] = child[3].xinterval;
      intervalB[3] = child[3].yinterval;
      coordA = particles.x;
      coordB = particles.y;
      startpointA = tp_x0;
      endpointA = tp_x1;
      startpointB = tp_y0;
      endpointB = tp_y1;
      break;

    case 3: // XYZ = 011, A is Y and B is Z
      yL *= 0.5;
      zL *= 0.5;
      intervalA[0] = child[0].yinterval;
      intervalB[0] = child[0].zinterval;
      intervalA[1] = child[1].yinterval;
      intervalB[1] = child[1].zinterval;
      intervalA[2] = child[2].yinterval;
      intervalB[2] = child[2].zinterval;
      intervalA[3] = child[3].yinterval;
      intervalB[3] = child[3].zinterval;
      coordA = particles.y;
      coordB = particles.z;
      startpointA = tp_y0;
      endpointA = tp_y1;
      startpointB = tp_z0;
      endpointB = tp_z1;
      break;

    case 5: // XYZ = 101, A is Z and B is X
      zL *= 0.5;
      xL *= 0.5;
      intervalA[0] = child[0].zinterval;
      intervalB[0] = child[0].xinterval;
      intervalA[1] = child[1].zinterval;
      intervalB[1] = child[1].xinterval;
      intervalA[2] = child[2].zinterval;
      intervalB[2] = child[2].xinterval;
      intervalA[3] = child[3].zinterval;
      intervalB[3] = child[3].xinterval;
      coordA = particles.z;
      coordB = particles.x;
      startpointA = tp_z0;
      endpointA = tp_z1;
      startpointB = tp_x0;
      endpointB = tp_x1;
      break;

    default:
      break;
  }

  midpointA = 0.5 * (startpointA + endpointA);
  midpointB = 0.5 * (startpointB + endpointB);

  // Child 0 ends with mid point on axis A, and ends with mid point on axis B
  intervalA[0][1] = midpointA;
  intervalB[0][1] = midpointB;

  // Child 1 begins with mid point on axis A, and ends with mid point on axis B
  intervalA[1][0] = midpointA;
  intervalB[1][1] = midpointB;

  // Child 2 ends with mid point on axis A, and begins with mid point on axis B
  intervalA[2][1] = midpointA;
  intervalB[2][0] = midpointB;

  // Child 3 begins with mid point on axis A, and begins with mid point on axis B
  intervalA[3][0] = midpointA;
  intervalB[3][0] = midpointB;

  double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
  double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

  for (int i = 0; i < 4; i++) {
    child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
    child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
    child[i].zc = 0.5 * (child[i].zinterval[0] + child[i].zinterval[1]);
    child[i].MAC = MAC;
  }

  vector<size_t> v[4];
  size_t start = tree[panel_index].members[0];
  size_t end = tree[panel_index].members[1];
  size_t* addr_table = new size_t[end - start + 1];

  size_t index;
  for (index = start; index <= end; index++) {
    particles.index[index] = index;
    addr_table[index - start] = index;

    if (coordA[index] <= midpointA && coordB[index] <= midpointB)
      v[0].push_back(index);
    else if (coordA[index] > midpointA && coordB[index] <= midpointB)
      v[1].push_back(index);
    else if (coordA[index] <= midpointA && coordB[index] > midpointB)
      v[2].push_back(index);
    else if (coordA[index] > midpointA && coordB[index] > midpointB)
      v[3].push_back(index);
  }

  size_t seq = start;
  for (size_t j = 0; j < 4; j++) {
    size_t size = v[j].size();

    if (size >= 1) {
      for (size_t k = 0; k < size; k++) {
        if (k == 0)
          child[j].members[0] = seq;
        if (k == size - 1)
          child[j].members[1] = seq;

        index = v[j][k];
        /*
        // This is very slow
        size_t pos;
        for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++) {
          if (particles.index[pos] == index)
            break;
        }
        Swap(pos, seq);
        */
        // This uses an address table
        size_t pos = addr_table[index - start];
        size_t out = particles.index[seq];
        Swap(pos, seq);
        addr_table[index - start] = seq;
        addr_table[out - start] = pos;

        seq++;
      }

      node_count++;
      tree[panel_index].children.push_back(node_count - 1);
      tree.push_back(child[j]);
      v[j].clear();
    }
  }

  delete[] addr_table;
}

void split_8(size_t panel_index)
{
  panel child[8];

/*
      ^ axis y
      |
  end -------------------
      |        |        |
      |Child 2 |Child 3 |
  mid |------------------
      |        |        |
start |Child 0 |Child 1 |
-----------------------------> axis x (lower z level)
    start     mid      end

      ^ axis y
      |
  end -------------------
      |        |        |
      |Child 6 |Child 7 |
  mid |------------------
      |        |        |
start |Child 4 |Child 5 |
-----------------------------> axis x (upper z level)
    start     mid      end
*/

  double tp_x0 = tree[panel_index].xinterval[0];
  double tp_x1 = tree[panel_index].xinterval[1];
  double tp_y0 = tree[panel_index].yinterval[0];
  double tp_y1 = tree[panel_index].yinterval[1];
  double tp_z0 = tree[panel_index].zinterval[0];
  double tp_z1 = tree[panel_index].zinterval[1];

  double xL = 0.5 * (tp_x1 - tp_x0);
  double yL = 0.5 * (tp_y1 - tp_y0);
  double zL = 0.5 * (tp_z1 - tp_z0);

  double midpointx = 0.5 * (tp_x0 + tp_x1);
  double midpointy = 0.5 * (tp_y0 + tp_y1);
  double midpointz = 0.5 * (tp_z0 + tp_z1);

  child[0].xinterval[0] = tp_x0;
  child[0].xinterval[1] = midpointx;
  child[0].yinterval[0] = tp_y0;
  child[0].yinterval[1] = midpointy;
  child[0].zinterval[0] = tp_z0;
  child[0].zinterval[1] = midpointz;

  child[1].xinterval[0] = midpointx;
  child[1].xinterval[1] = tp_x1;
  child[1].yinterval[0] = tp_y0;
  child[1].yinterval[1] = midpointy;
  child[1].zinterval[0] = tp_z0;
  child[1].zinterval[1] = midpointz;

  child[2].xinterval[0] = tp_x0;
  child[2].xinterval[1] = midpointx;
  child[2].yinterval[0] = midpointy;
  child[2].yinterval[1] = tp_y1;
  child[2].zinterval[0] = tp_z0;
  child[2].zinterval[1] = midpointz;

  child[3].xinterval[0] = midpointx;
  child[3].xinterval[1] = tp_x1;
  child[3].yinterval[0] = midpointy;
  child[3].yinterval[1] = tp_y1;
  child[3].zinterval[0] = tp_z0;
  child[3].zinterval[1] = midpointz;

  child[4].xinterval[0] = tp_x0;
  child[4].xinterval[1] = midpointx;
  child[4].yinterval[0] = tp_y0;
  child[4].yinterval[1] = midpointy;
  child[4].zinterval[0] = midpointz;
  child[4].zinterval[1] = tp_z1;

  child[5].xinterval[0] = midpointx;
  child[5].xinterval[1] = tp_x1;
  child[5].yinterval[0] = tp_y0;
  child[5].yinterval[1] = midpointy;
  child[5].zinterval[0] = midpointz;
  child[5].zinterval[1] = tp_z1;

  child[6].xinterval[0] = tp_x0;
  child[6].xinterval[1] = midpointx;
  child[6].yinterval[0] = midpointy;
  child[6].yinterval[1] = tp_y1;
  child[6].zinterval[0] = midpointz;
  child[6].zinterval[1] = tp_z1;

  child[7].xinterval[0] = midpointx;
  child[7].xinterval[1] = tp_x1;
  child[7].yinterval[0] = midpointy;
  child[7].yinterval[1] = tp_y1;
  child[7].zinterval[0] = midpointz;
  child[7].zinterval[1] = tp_z1;

  double sq_r = 0.25 * (xL * xL + yL * yL + zL * zL); // r^2
  double MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

  for (int i = 0; i < 8; i++) {
    child[i].xc = 0.5 * (child[i].xinterval[0] + child[i].xinterval[1]);
    child[i].yc = 0.5 * (child[i].yinterval[0] + child[i].yinterval[1]);
    child[i].zc = 0.5 * (child[i].zinterval[0] + child[i].zinterval[1]);
    child[i].MAC = MAC;
  }

  vector<size_t> v[8];
  size_t start = tree[panel_index].members[0];
  size_t end = tree[panel_index].members[1];
  size_t* addr_table = new size_t[end - start + 1];

  size_t index;
  for (index = start; index <= end; index++) {
    particles.index[index] = index;
    addr_table[index - start] = index;

    if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
        particles.z[index] <= midpointz)
      v[0].push_back(index);
    else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
             particles.z[index] <= midpointz)
      v[1].push_back(index);
    else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
             particles.z[index] <= midpointz)
      v[2].push_back(index);
    else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
             particles.z[index] <= midpointz)
      v[3].push_back(index);
    else if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
             particles.z[index] > midpointz)
      v[4].push_back(index);
    else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
             particles.z[index] > midpointz)
      v[5].push_back(index);
    else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
             particles.z[index] > midpointz)
      v[6].push_back(index);
    else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
             particles.z[index] > midpointz)
      v[7].push_back(index);
  }

  size_t seq = start;
  for (size_t j = 0; j < 8; j++) {
    size_t size = v[j].size();

    if (size >= 1) {
      for (size_t k = 0; k < size; k++) {
        if (k == 0)
          child[j].members[0] = seq;
        if (k == size - 1)
          child[j].members[1] = seq;

        index = v[j][k];
        /*
        // This is very slow
        size_t pos;
        for (pos = tree[panel_index].members[0]; pos <= tree[panel_index].members[1]; pos++) {
          if (particles.index[pos] == index)
            break;
        }
        Swap(pos, seq);
        */
        // This uses an address table
        size_t pos = addr_table[index - start];
        size_t out = particles.index[seq];
        Swap(pos, seq);
        addr_table[index - start] = seq;
        addr_table[out - start] = pos;

        seq++;
      }

      node_count++;
      tree[panel_index].children.push_back(node_count - 1);
      tree.push_back(child[j]);
      v[j].clear();
    }
  }

  delete[] addr_table;
}

void split_tree_node(size_t panel_index)
{
  double xL = tree[panel_index].xinterval[1] - tree[panel_index].xinterval[0];
  double yL = tree[panel_index].yinterval[1] - tree[panel_index].yinterval[0];
  double zL = tree[panel_index].zinterval[1] - tree[panel_index].zinterval[0];
  double L_max = xL;

  if (yL > L_max)
    L_max = yL;

  if (zL > L_max)
    L_max = zL;

  int XYZ_flag = 0;
  const double ratio = sqrt(2.0);

  if (xL * ratio > L_max)
    XYZ_flag += 4;

  if (yL * ratio > L_max)
    XYZ_flag += 2;

  if (zL * ratio > L_max)
    XYZ_flag += 1;

  switch (XYZ_flag) {
    case 1: // XYZ = 001, split along Z
    case 2: // XYZ = 010, split along Y
    case 4: // XYZ = 100, split along X
#ifdef TREE_STAT
      split2_times++;
#endif
      split_2(panel_index, XYZ_flag);
      break;

    case 3: // XYZ = 011, split along Y and Z
    case 5: // XYZ = 101, split along Z and X
    case 6: // XYZ = 110, split along X and Y
#ifdef TREE_STAT
      split4_times++;
#endif
      split_4(panel_index, XYZ_flag);
      break;

    case 7: // XYZ = 111, split along X, Y, and Z
#ifdef TREE_STAT
      split8_times++;
#endif
      split_8(panel_index);
      break;

    default:
      break;
  }
}

void build_tree_3D_Recursive(size_t panel_index, int level)
{
#ifdef TREE_STAT
  if (level > max_level)
    max_level = level;

  level_cnt[level]++;
#endif

  size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;
  if (n >= (size_t)N0) {
    split_tree_node(panel_index);

    for (size_t i = 0; i < tree[panel_index].children.size(); i++) {
      size_t panel_index_new = tree[panel_index].children[i];
      build_tree_3D_Recursive(panel_index_new, level + 1);
    }
  }
  else {
    leaf.push_back(panel_index);

#ifdef TREE_STAT
    double tp_x0 = tree[panel_index].xinterval[0];
    double tp_x1 = tree[panel_index].xinterval[1];
    double tp_y0 = tree[panel_index].yinterval[0];
    double tp_y1 = tree[panel_index].yinterval[1];
    double tp_z0 = tree[panel_index].zinterval[0];
    double tp_z1 = tree[panel_index].zinterval[1];

    double xL = tp_x1 - tp_x0;
    double yL = tp_y1 - tp_y0;
    double zL = tp_z1 - tp_z0;

    double diameter = sqrt(xL * xL + yL * yL + zL * zL);
    if (diameter > max_leaf_diameter)
      max_leaf_diameter = diameter;
    if (diameter < min_leaf_diameter)
      min_leaf_diameter = diameter;

    int pars = tree[panel_index].members[1] - tree[panel_index].members[0] + 1;
    if (pars > max_leaf_pars)
      max_leaf_pars = pars;
    if (pars < min_leaf_pars)
      min_leaf_pars = pars;

    double Lmax = xL;
    if (yL > Lmax)
      Lmax = yL;
    if (zL > Lmax)
      Lmax = zL;

    double Lmin = xL;
    if (yL < Lmin)
      Lmin = yL;
    if (zL < Lmin)
      Lmin = zL;

    double leaf_ratio = Lmax / Lmin;
    if (leaf_ratio > max_leaf_ratio)
      max_leaf_ratio = leaf_ratio;
    if (leaf_ratio < min_leaf_ratio)
      min_leaf_ratio = leaf_ratio;
#endif
  }
}

void Panel_Moment_B(size_t panel_index, double wx1[P + 1], double wy1[P + 1],
                    double wz1[P + 1], double m[][Pflat])
{
  double t1[P + 1];
  double t2[P + 1];
  double t3[P + 1];
  double wx0[P + 1];
  double wy0[P + 1];
  double wz0[P + 1];

  int i, j, k, kk;
  for (i = 0; i < P + 1; i++) {
    t1[i] = tree[panel_index].t1[i];
    t2[i] = tree[panel_index].t2[i];
    t3[i] = tree[panel_index].t3[i];
    wx0[i] = tree[panel_index].wx0[i];
    wy0[i] = tree[panel_index].wy0[i];
    wz0[i] = tree[panel_index].wz0[i];
  }

  size_t tp0 = tree[panel_index].members[0];
  size_t tp1 = tree[panel_index].members[1];
  size_t tp_j;

  double tempx[P + 1];
  double tempy[P + 1];
  double tempz[P + 1];
  double tempx1[P + 1];
  double tempy1[P + 1];
  double tempz1[P + 1];
  int a1exactind, a2exactind, a3exactind;

  double x, y, z;
  double dx, dy, dz;
  double sum1;
  double sum2;
  double sum3;

  double D, lam;
  double temp11, temp21, temp12, temp22;
  int index;

  kk = -1;
  for (i = 0; i < P + 1; i++) {
    for (j = 0; j < P + 1; j++) {
      for (k = 0; k < P + 1; k++) {
        kk++;
        for (index = 0; index < 8; index++)
          m[index][kk] = 0;
      }
    }
  }

  for (tp_j = tp0; tp_j <= tp1; tp_j++) {
    x = particles.x[tp_j];
    y = particles.y[tp_j];
    z = particles.z[tp_j];
    lam = lambda[tp_j];

    a1exactind = -1;
    a2exactind = -1;
    a3exactind = -1;

    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;

    for (j = 0; j < P + 1; j++) {
      dx = x - t1[j];
      dy = y - t2[j];
      dz = z - t3[j];

      if (fabs(dx) <= DBL_MIN)
        a1exactind = j;
      else {
        tempx[j] = wx0[j] / dx + wx1[j] / (dx * dx);
        tempx1[j] = wx1[j] / dx;
        sum1 += tempx[j];
      }

      if (fabs(dy) <= DBL_MIN)
        a2exactind = j;
      else {
        tempy[j] = wy0[j] / dy + wy1[j] / (dy * dy);
        tempy1[j] = wy1[j] / dy;
        sum2 += tempy[j];
      }

      if (fabs(dz) <= DBL_MIN)
        a3exactind = j;
      else {
        tempz[j] = wz0[j] / dz + wz1[j] / (dz * dz);
        tempz1[j] = wz1[j] / dz;
        sum3 += tempz[j];
      }
    }

    if (a1exactind > -1) {
      sum1 = 1.0;

      for (j = 0; j < P + 1; j++) {
        tempx[j] = 0.0;
        tempx1[j] = 0.0;
      }

      tempx[a1exactind] = 1.0;
      tempx1[a1exactind] = 0;
    }

    if (a2exactind > -1) {
      sum2 = 1.0;

      for (j = 0; j < P + 1; j++) {
        tempy[j] = 0.0;
        tempy1[j] = 0.0;
      }

      tempy[a2exactind] = 1.0;
      tempy1[a2exactind] = 0;
    }

    if (a3exactind > -1) {
      sum3 = 1.0;

      for (j = 0; j < P + 1; j++) {
        tempz[j] = 0.0;
        tempz1[j] = 0.0;
      }

      tempz[a3exactind] = 1.0;
      tempz1[a3exactind] = 0;
    }

    D = 1.0 / (sum1 * sum2 * sum3) * lam;
    kk = -1;

    for (i = 0; i < P + 1; i++) {
      for (j = 0; j < P + 1; j++) {
        for (k = 0; k < P + 1; k++) {
          kk++;
          temp11 = tempx[i] * tempy[j] * D;
          temp21 = tempx1[i] * tempy[j] * D;
          temp12 = tempx[i] * tempy1[j] * D;
          temp22 = tempx1[i] * tempy1[j] * D;
          m[0][kk] += temp11 * tempz[k];
          m[1][kk] += temp21 * tempz[k];
          m[2][kk] += temp12 * tempz[k];
          m[3][kk] += temp11 * tempz1[k];
          m[4][kk] += temp22 * tempz[k];
          m[5][kk] += temp12 * tempz1[k];
          m[6][kk] += temp21 * tempz1[k];
          m[7][kk] += temp22 * tempz1[k];
        }
      }
    }
  }
}

double Call_Treecode(double x, double y, double z, int panel_index)
{
  double velocity = 0.0;
  double a, ap1, ap2, ap3, ap12, ap23, ap13, ap123;
  double d1, d2, d3;
  double R2, R, Rinv, Rinv3, Rinv5, Rinv7;

  double t1[P + 1];
  double t2[P + 1];
  double t3[P + 1];

  int i, j, k, kk;
  for (i = 0; i < P + 1; i++) {
    t1[i] = tree[panel_index].t1[i];
    t2[i] = tree[panel_index].t2[i];
    t3[i] = tree[panel_index].t3[i];
  }

  kk = -1;
  double sq_d1, sq_d2, sq_d3;
  for (i = 0; i < P + 1; i++) {
    d1 = x - t1[i];
    sq_d1 = d1 * d1;
    for (j = 0; j < P + 1; j++) {
      d2 = y - t2[j];
      sq_d2 = d2 * d2;
      for (k = 0; k < P + 1; k++) {
        kk++;
        d3 = z - t3[k];
        sq_d3 = d3 * d3;
        R2 = sq_d1 + sq_d2 + sq_d3;
        R = sqrt(R2);
        Rinv = 1.0 / R;
        Rinv3 = Rinv * Rinv * Rinv;
        Rinv5 = Rinv3 * Rinv * Rinv;
        Rinv7 = Rinv5 * Rinv * Rinv;

        a = Rinv;
        ap1 = d1 * Rinv3;
        ap2 = d2 * Rinv3;
        ap3 = d3 * Rinv3;
        ap12 = 3.0 * d1 * d2 * Rinv5;
        ap23 = 3.0 * d2 * d3 * Rinv5;
        ap13 = 3.0 * d1 * d3 * Rinv5;
        ap123 = 15.0 * d1 * d2 * d3 * Rinv7;

        velocity += a * tree[panel_index].moments[0][kk] + ap1* tree[panel_index].moments[1][kk] +
                    ap2 * tree[panel_index].moments[2][kk] + ap3 * tree[panel_index].moments[3][kk] +
                    ap12 * tree[panel_index].moments[4][kk] + ap23 * tree[panel_index].moments[5][kk] +
                    ap13 * tree[panel_index].moments[6][kk] + ap123 * tree[panel_index].moments[7][kk];
      }
    }
  }

  return velocity;
}

double Call_Ds(size_t limit_1, size_t limit_2, size_t particle_index, double p_x, double p_y, double p_z)
{
  double xx, yy, zz;
  double R2, R, Rinv;

  double velocity = 0.0;

  for (size_t jj = limit_1; jj <= limit_2; jj++) {
    xx = p_x - particles.x[jj];
    yy = p_y - particles.y[jj];
    zz = p_z - particles.z[jj];

    R2 = xx * xx + yy * yy + zz * zz;
    R = sqrt(R2);
    if (jj != particle_index) {
      Rinv = 1.0 / R;
      velocity += Rinv * lambda[jj];
    }
  }

  return velocity;
}

double Compute_RBF(size_t particle_index, size_t panel_index)
{
  // Input :
  //         particle_index
  //         panel_index
  // Output :
  //         velocity in 1D
  double velocity = 0.0;

  double p_x = 0.0;
  double p_y = 0.0;
  double p_z = 0.0;
  double xc = 0.0;
  double yc = 0.0;
  double zc = 0.0;
  size_t limit_1;
  size_t limit_2;
  double R_sq;

  limit_1 = tree[panel_index].members[0];
  limit_2 = tree[panel_index].members[1];

  p_x = particles.x[particle_index];
  p_y = particles.y[particle_index];
  p_z = particles.z[particle_index];

  xc = tree[panel_index].xc;
  yc = tree[panel_index].yc;
  zc = tree[panel_index].zc;

  double tpx = p_x - xc;
  double tpy = p_y - yc;
  double tpz = p_z - zc;

  R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

  if (tree[panel_index].MAC < R_sq)
    velocity = Call_Treecode(p_x, p_y, p_z, panel_index);
  else {
    if (limit_2 - limit_1 < N0) // Otherwise, if cluster is a leaf, use direct sum
      velocity = Call_Ds(limit_1, limit_2, particle_index, p_x, p_y, p_z);
    else { // Otherwise, if cluster is not a leaf, look at children
      size_t length = tree[panel_index].children.size();
      for (size_t i = 0; i < length; i++) {
        size_t index = tree[panel_index].children[i];
        double temp_result = Compute_RBF(particle_index, index);
        velocity += temp_result;
      }
    }
  }

  return velocity;
}

//*****************************************************************************//
void Chev_Points(double t[])
{
  double h = pi / P;
  for (int i = 0; i < P + 1; i++)
    t[i] = cos(i * h); // Chebyshev interpolation points [-1, 1]
}

void Cluster_Chev_PointsWeights(size_t tree_size, double t[], double wx1[P + 1],
                                double wy1[P + 1], double wz1[P + 1])
{
  double x1, x2, y1, y2, z1, z2;
  size_t tree_index;
  double wx0[P + 1] = {0.0};
  double wy0[P + 1] = {0.0};
  double wz0[P + 1] = {0.0};

  wx0[0] = P * P / 3.0 + 1.0 / 6.0;
  wx0[P] = -wx0[0];
  wy0[0] = wx0[0];
  wy0[P] = -wy0[0];
  wz0[0] = wx0[0];
  wz0[P] = -wz0[0];

  double s;
  for (int j = 1; j < P; j++) {
    s = j * pi / P;
    wx0[j] = -cos(s) / (2.0 * sin(s) * sin(s));
    wy0[j] = wx0[j];
    wz0[j] = wx0[j];
  }

  double dx, dy, dz;
  for (tree_index = 0; tree_index < tree_size; tree_index++) {
    x1 = tree[tree_index].xinterval[0];
    x2 = tree[tree_index].xinterval[1];
    y1 = tree[tree_index].yinterval[0];
    y2 = tree[tree_index].yinterval[1];
    z1 = tree[tree_index].zinterval[0];
    z2 = tree[tree_index].zinterval[1];
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    for (int i = 0; i < P + 1; i++) { // Map to the cluster
      tree[tree_index].t1[i] = x1 + (t[i] + 1.0) * 0.5 * dx;
      tree[tree_index].t2[i] = y1 + (t[i] + 1.0) * 0.5 * dy;
      tree[tree_index].t3[i] = z1 + (t[i] + 1.0) * 0.5 * dz;

      tree[tree_index].wx0[i] = -2.0 * wx1[i] * 2.0 * wx0[i] /dx;
      tree[tree_index].wy0[i] = -2.0 * wy1[i] * 2.0 * wy0[i] /dy;
      tree[tree_index].wz0[i] = -2.0 * wz1[i] * 2.0 * wz0[i] /dz;
    }
  }
}

int main()
{
  cout << " ===== Hermite 1/r ===========" << endl;
  cout << "P is " << P << endl;
  cout << "numpars is " << numpars << endl;
  cout << "theta is " << sqrt(sq_theta) << endl;
  cout << "N0 is " << N0 << endl;

  char point_data_file[64] = {0};
  sprintf(point_data_file, "./rand_%d.txt", numpars);

  FILE *fp = fopen(point_data_file, "r");
  if (fp == NULL) {
    cerr << "Cannot open point data file " << point_data_file << "!" << endl;
    return 1;
  }

  int counter = 0;
  double x1, x2, x3;
  while (fscanf(fp, "%lf %lf %lf", &x1, &x2, &x3) == 3) {
    particles.x[counter] = x1;
    particles.y[counter] = x2;
    particles.z[counter] = x3;
    particles.index[counter] = -1;
    particles.old_index[counter] = counter;
    counter++;
    if (counter == numpars)
      break;
  }
  fclose(fp);

  if (counter < numpars) {
    cerr << "Less lines of point data were read! counter = " << counter << endl;
    return 1;
  }

  // Find bounds of Cartesian box enclosing the particles
  xyzminmax[0] = minval(particles.x);
  xyzminmax[1] = maxval(particles.x);
  xyzminmax[2] = minval(particles.y);
  xyzminmax[3] = maxval(particles.y);
  xyzminmax[4] = minval(particles.z);
  xyzminmax[5] = maxval(particles.z);

  cout << endl;
  cout << "xmin = " << xyzminmax[0] << ", xmax = " << xyzminmax[1] << endl;
  cout << "ymin = " << xyzminmax[2] << ", ymax = " << xyzminmax[3] << endl;
  cout << "zmin = " << xyzminmax[4] << ", zmax = " << xyzminmax[5] << endl;
  cout << endl;

  lambda = new double[numpars];

  char lambda_data_file[64] = {0};
  sprintf(lambda_data_file, "./lambda_%d.txt", numpars);

  fp = fopen(lambda_data_file, "r");
  if (fp == NULL) {
    cerr << "Cannot open lambda data file " << lambda_data_file << "!" << endl;
    return 1;
  }

  counter = 0;
  while (fscanf(fp, "%lf", &x1) == 1) {
    lambda[counter] = x1;
    counter++;
    if (counter == numpars)
      break;
  }
  fclose(fp);

  if (counter < numpars) {
    cerr << "Less lines of lambda data were read! counter = " << counter << endl;
    return 1;
  }

  velo = new double[numpars];
  memset(velo, 0, numpars * sizeof(double));

  velo_true = new double[numpars];
  memset(velo_true, 0, numpars * sizeof(double));

  //***************** Set up tree *******************************
  long total_time, build_tree_time, moment_cpu_time,treecode_cpu_time, ds_cpu_time;
  long build_tree_end_time, moment_cpu_end_time;
    
  total_time = getTickCount();

  build_tree_init();
  build_tree_3D_Recursive(0, 0);

  build_tree_end_time = getTickCount();

  //***************** Compute moment for each panel **************
    
  moment_cpu_time = getTickCount();
  double t[P + 1];
  Chev_Points(t); // Cheby points at [-1, 1]

  double wx1[P + 1] = {0.0};
  double wy1[P + 1] = {0.0};
  double wz1[P + 1] = {0.0};

  wx1[0] = 0.25;
  wx1[P] = 0.25;
  wy1[0] = 0.25;
  wy1[P] = 0.25;
  wz1[0] = 0.25;
  wz1[P] = 0.25;

  for (int j = 1; j < P; j++) {
    wx1[j] = 1.0;
    wy1[j] = 1.0;
    wz1[j] = 1.0;
  }

  size_t tree_size = tree.size();
  Cluster_Chev_PointsWeights(tree_size, t, wx1, wy1, wz1);

  // Skip root
  for (size_t i = 1; i < tree_size; i++)
    Panel_Moment_B(i, wx1, wy1, wz1, tree[i].moments);

  moment_cpu_end_time = getTickCount();
  //***************** Compute Velocity ***************************
  for (int particle_index = 0; particle_index < numpars; particle_index++)
    velo[particle_index] = Compute_RBF(particle_index, 0);

  treecode_cpu_time = getTickCount() - total_time;

  cout << "build_tree_time " << build_tree_end_time - total_time << endl;
  cout << "moment_time " << moment_cpu_end_time - moment_cpu_time << endl;
  cout << "treecode_cpu_time " << treecode_cpu_time << endl;

#ifdef TREE_STAT
  cout << endl;
  cout << "tree size = " << tree.size() << ", leaf count = " << leaf.size() << endl;
  cout << "max level = " << max_level << endl;
 
  for (int i = 0; i <= max_level; i++)
    cout << "level " << i << ": " << level_cnt[i] << " panel(s)" << endl;

  cout << "split2_times = " << split2_times << ", split4_times = " << split4_times << ", split8_times = " << split8_times << endl;
  cout << "max_leaf_pars = " << max_leaf_pars << ", min_leaf_pars = " << min_leaf_pars << endl;
  cout << "max_leaf_ratio = " << max_leaf_ratio << ", min_leaf_ratio = " << min_leaf_ratio << endl;
  cout << "max_leaf_diameter = " << max_leaf_diameter << ", min_leaf_diameter = " << min_leaf_diameter << endl;
  cout << endl;
#endif

  //***************** Compute true Velocity with Direct summation ***********
  ds_cpu_time = getTickCount(); 

  compute_Ds();

  ds_cpu_time = getTickCount() - ds_cpu_time;

  cout << "ds_cpu_time " << ds_cpu_time << endl;

  cout << "\nDisplay first 5 potential values calculated by compute_Ds" << endl;
  for (int i = 0; i < 5; i++)
    cout << setprecision(25) << "velo_true[" << i << "] = " << velo_true[i] << " " << endl;

  cout << "\nDisplay first 5 potential values calculated by treecode" << endl;
  for (int i = 0; i < 5; i++)
    cout << setprecision(25) << "velo[" << i << "] = " << velo[i] << " " << endl;

  //******** Error: L_infty *************************
  double E = 0.0;
  double temp_n = 0.0;
  double temp_d = 0.0;
  double max_n = 0.0;
  double max_d = 0.0;

  for (int i = 0; i < numpars; i++) {
    temp_n = (velo[i] - velo_true[i]) * (velo[i] - velo_true[i]);
    temp_d = velo_true[i] * velo_true[i];

    temp_n = sqrt(temp_n);
    temp_d = sqrt(temp_d);
    if (temp_n > max_n)
      max_n = temp_n;
    if (temp_d > max_d)
      max_d = temp_d;
    temp_d = 0.0;
    temp_n = 0.0;
  }

  E = max_n / max_d;

  cout << endl;
  cout << setprecision(25) << "E is " << E << endl;
  cout << setprecision(25) << "max_n is " << max_n << endl;
  cout << setprecision(25) << "max_d is " << max_d << endl;

  //******** Error: extend from RBF paper *******************
  double err2_ex = 0.0;
  double sum_d_ex = 0.0;
  double sum_n_ex = 0.0;
  for (int i = 0; i < numpars; i++) {
    sum_n_ex += (velo[i] - velo_true[i]) * (velo[i] - velo_true[i]);
    sum_d_ex += velo_true[i] * velo_true[i];
  }

  err2_ex = sqrt(sum_n_ex / sum_d_ex);

  cout << endl;
  cout << setprecision(25) << "E_2_ex is " << err2_ex << endl;

  //*********************************************************
  delete[] lambda;
  delete[] velo;
  delete[] velo_true;

  cout << "\nDone" << endl;

  return 0;
}
