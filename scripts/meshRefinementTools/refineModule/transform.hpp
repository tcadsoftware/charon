
namespace lcm_lib
{

  class Transform
    {

    public:

      Transform();

      /*
      Transform(double,double,double,double,double,double,
		double,double,double);
      */

      ~Transform();

      void transform(double ,double ,double ,double ,
		     double ,double ,double ,double ,
		     double ,double& ,double& ,double& ,
		     double& ,double& ,double& );

      void planar_coords(double&, double&, double, double, double);
      void deplanar_coords(double, double, double&, double&, double&);


    private:

      int xform_check;
      double dir_cos[3][3];
      double translation[3];


      void normalize(double[], int);

      void cross_product(double[],double[],double[], int);

      double dot_product(double[], double[], int);

      void get_dir_cos(double[],double[],double[],double[],double[],
		       double[], int);

      void mat_vec(double[], double[], int, int);
      void trans_mat_vec(double[], double[], int, int);


    };

  //Some other functions I want to keep in my namespace

  double tri_area(double,double,double,double,double,double);

  bool tri_check(double,double,double,double,
		 double,double,double,double);

}
