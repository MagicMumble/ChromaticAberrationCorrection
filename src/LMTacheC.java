import static java.lang.Math.*;

public class LMTacheC {
    private boolean clr;
    private int wi, he;
    private int xbegin, xend, ybegin, yend;
    private Image im;
    static final double EPSILON_KERNEL = 1E-9;

    public LMTacheC(Image Im, double cx, double cy, double delta, boolean tache_color, int img_width, int img_height) {
        clr = tache_color;
        wi = img_width;
        he = img_height;
        xbegin = (int) (cx + 0.5 - delta);
        xend = (int) (cx + 0.5 + delta);
        ybegin = (int) (cy + 0.5 - delta);
        yend = (int) (cy + 0.5 + delta);
        im = Im;
    }

    public double[] minus(double[] mas1, double[] mas2){
        double[] mas = new double[mas1.length];
        for (int i=0; i< mas1.length; i++)
            mas[i]=mas1[i]-mas2[i];
        return mas;
    }

    public double[] plus(double[] mas1, double[] mas2){
        double[] mas = new double[mas1.length];
        for (int i=0; i< mas1.length; i++)
            mas[i]=mas1[i]+mas2[i];
        return mas;
    }

    public double qnorm(double[] mas) {
        double q = 0.0;
        for (int i = mas.length; i >= 0; i--)
            q += mas[i]*mas[i];
        return q;
    }

    public void modelData(double[] P, double[] ymodel) {
        int idx = 0;
        for (int i = ybegin; i <= yend; i++) {
            for (int j = xbegin; j <= xend; j++) {
                ymodel[idx] = model_luminance_alternative(j, i, wi, he, P, clr);
                idx++;
            }
        }
    }

    public double model_luminance_alternative(int x, int y, int wi, int he, double[] P, boolean clr)
    {
        double err = 0, lum = 0;
        double dist = -512;
        double l1=P[0], l2=P[1], theta=P[2], tu=P[3], tv=P[4], r1=P[5], p=P[6], r2=P[7], vh=P[8], vb=P[9], c=P[10];
        if (x >= 1 && x <= wi-2 && y >= 1 && y <= he-2)
        {
            double coordX = l1 * (x-tu)* cos(theta) - l1 * (y-tv)* sin(theta);
            double coordY = l2 * (x-tu)* sin(theta) + l2 * (y-tv)* cos(theta);
            double d = coordX*coordX + coordY*coordY;
            if (d < 2.0*2.0)
                dist = Math.sqrt(d);
        }

        double d = dist;
        if (d == -512) // does not belong to circle - return the 'not belongin color'.
        {
            lum = 255;
            vh = 255;
            vb = 0;
            if (!clr)
                lum = vb + lum;
            else
                lum = vh - lum;
            return lum;
        }
        if (!clr)
            d = -d+c;
        else
            d = d-c;
        r1 = 0;//abs(r1);
        r2 = 0; //abs(r2);
        double y20 = -1.0;
        double x21 = y20/p;
        double x20 = x21;
        if (d >= x20)
            lum = 0;
        else
        {
            double y10 = 1.0;
            double x11 = y10/p;
            double x10 = x11;
            if (d <= x10)
                lum = 1;
            else
            {
                lum = p*d;
                lum = 0.5*lum+0.5;
            }
        }

        lum = lum * (P[8] - P[9]);
        if (!clr)
            lum = P[9] + lum;
        else
            lum = P[8] - lum;
        return lum;
    }

    void jacobian_alternative(Matrix J, double[] P, int wi, int he, boolean clr, int xbegin, int xend, int ybegin, int yend)
    {
        double l1=P[0], l2=P[1], theta=P[2], tu=P[3], tv=P[4], r1=P[5], p=P[6], r2=P[7], vh=P[8], vb=P[9], c=P[10];
        int idx = 0;
        for (int v = ybegin; v <= yend; v++) {
            for (int u = xbegin; u <= xend; u++) {

                for (int k = 0; k < P.length; k++) // initialize
                    J.set(idx, k, 0);

                double x = 0, y = 0, d = -512;
                if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2) {
                    x = l1 * (u-tu)* cos(theta) - l1 * (v-tv)* sin(theta);
                    y = l2 * (u-tu)* sin(theta) + l2 * (v-tv)* cos(theta);
                    double dist = x*x + y*y;
                    if (dist <= 2.0*2.0)
                        d = Math.sqrt(dist);
                }

                if (d != -512) {
                    if (!clr) d = -d+c;
                    else d = d-c;
                    r1 = 0;//abs(r1);
                    r2 = 0;//abs(r2);

                    double y20 = -1.0;
                    double x21 = y20/p;
                    double x20 = x21;

                    double y10 = 1.0;
                    double x11 = y10/p;
                    double x10 = x11;

                    double s2=0;
                    if(d>=x20)
                        s2 = 0.0;
                    else if(d<=x10)
                        s2 = 1.0;
                    else //if (d<=x21)
                    {
                        s2=p*d;
                        s2 = 0.5*s2+0.5;}
                    if (d < x20 && d > x10) {
                        double dXdL1 = (u-tu)* cos(theta) - (v-tv)* sin(theta); // "-" - for black taches, because d = -d+1!
                        double dYdL2 = (u-tu)* sin(theta) + (v-tv)* cos(theta);
                        double dDdL1 = -x*dXdL1 / Math.sqrt(x*x + y*y);
                        double dDdL2 = -y*dYdL2 / Math.sqrt(x*x + y*y);
                        double dXdTheta = -l1*(sin(theta)*(u-tu) + cos(theta)*(v-tv));
                        double dYdTheta = l2*(cos(theta)*(u-tu) - sin(theta)*(v-tv));
                        double dDdTheta = -(x*dXdTheta + y*dYdTheta) / sqrt(x*x + y*y);
                        double dXdTu = -l1*cos(theta);
                        double dYdTu = -l2*sin(theta);
                        double dDdTu = -(x*dXdTu + y*dYdTu) / sqrt(x*x + y*y);
                        double dXdTv = l1*sin(theta);
                        double dYdTv = -l2*cos(theta);
                        double dDdTv = -(x*dXdTv + y*dYdTv) / sqrt(x*x + y*y);
                        double dDdC = 1;
                        J.set(idx, 0, p*dDdL1);
                        J.set(idx, 1, p*dDdL2);
                        J.set(idx, 2, p*dDdTheta);
                        J.set(idx, 3, p*dDdTu);
                        J.set(idx, 4, p*dDdTv);
                        //J(idx, 5) = 0;
                        J.set(idx, 6, d);
                        //J(idx, 7) = 0;
                        J.set(idx, 10, p*dDdC);
                        for (int k = 0; k < P.length; k++)
                            J.set(idx,k, J.get(idx,k)*0.5*(vh-vb));
                    }
                    J.set(idx, 8, s2);
                    J.set(idx, 9, 1-s2);
                }
                idx++;
            }
        }
    }

    public void modelJacobian(double[] P, Matrix J) {
        jacobian_alternative(J, P, wi, he, clr, xbegin, xend, ybegin, yend);
    }
}

