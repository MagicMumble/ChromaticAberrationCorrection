import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import static java.lang.Math.*;

public class Main {

    private static boolean color = false;
    private static double[] xR, yR, xGr, yGr, xB, yB, xGb, yGb, rR, rG, rB;

    public static void main(String[] args) throws IOException {
        if (args.length == 3) // only estimates and saves polynomial
            polyEstimation(args, color);
        else
            System.out.println("Enter .pgm file path and paths to two text files for writing down the polynomials");
    }

    public static void polyEstimation(String[] args, boolean color) throws IOException {
        System.out.println("\nPolynomial estimation... \n");
        String fnameRGB = args[0];
        String fnamePolyR = args[1];
        String fnamePolyB = args[2];
        int scale = 2;
        Image img_bayer = Image.read_pgm(fnameRGB);
        int wi = img_bayer.xsize, he = img_bayer.ysize;
        int wiRB = wi/2, heRB = he/2;
        int wiG = wiRB*scale, heG = heRB*scale;
        Image imgR = new Image(wiRB, heRB, 255);
        Image imgG = new Image(wiG, heG, 255);
        Image imgB = new Image(wiRB, heRB, 255);
        raw2rgb(img_bayer, imgR, imgG, imgB);
        keypnts_circle(imgR, imgG, imgB, scale, color);
        double[] paramsXR;
        double[] paramsYR;
        double[] paramsXB;
        double[] paramsYB;
        int degX = 11, degY = 11;
        double xp = (double)imgG.xsize/2+0.2, yp = (double)imgG.ysize/2+0.2;
        paramsXR=get_polynom1(xR, yR, xGr, yGr, degX, degY, xp, yp);
        paramsYR=get_polynom2(xR, yR, xGr, yGr, degX, degY, xp, yp);
        paramsXB=get_polynom1(xB, yB, xGb, yGb, degX, degY, xp, yp);
        paramsYB=get_polynom2(xB, yB, xGb, yGb, degX, degY, xp, yp);

        save_poly(fnamePolyR, paramsXR, paramsYR, degX, degY);
        save_poly(fnamePolyB, paramsXB, paramsYB, degX, degY);
    }

    public static void printMono(FileWriter pfile, double mono, int degX, int degY) throws IOException {
        if (degX != 0 && degY != 0) {
            if (mono > 0)	pfile.write("+ " + mono +" * x^" + degX + " * y^" + degY + "\n");
            else			pfile.write("- " + (-1)*mono +" * x^" + degX + " * y^" + degY + "\n"); }
        else if (degX == 0 && degY == 0) {
            if (mono > 0)	pfile.write("+ "+mono+"\n");
            else			pfile.write("- "+(-1)*mono+"\n"); }
        else if (degX == 0) {
            if (mono > 0)	pfile.write("+ "+mono + " * y^" + degY + "\n");
            else			pfile.write("- "+(-1)*mono + " * y^" + degY + "\n"); }
        else {
            if (mono > 0)	pfile.write("+ "+mono + " * x^" + degX + "\n");
            else			pfile.write("- "+(-1)*mono + " * x^" + degX + "\n"); }
    }

    public static void save_poly(String fname, double[] paramsX, double[] paramsY,  int degX, int degY) throws IOException {
        FileWriter writer = new FileWriter(fname, false);
        writer.write("# polyX(x,y): \n");
        int idx = 0;
        for (int i = degX; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                printMono(writer, paramsX[idx], i-j, j);
                idx++; 	} }
        idx = 0;
        writer.write("# polyY(x,y): \n");
        for (int i = degY; i >= 0; i--) {
            for (int j = 0; j <= i; j++) {
                printMono(writer, paramsY[idx], i-j, j);
                idx++; } }
        writer.close();
    }

    public static double[] get_polynom1(double[] xF, double[] yF, double[] xGf, double[] yGf,
                     int degX, int degY, double xp, double yp)
    {
        System.out.println("Obtaining correction polynomials... ");
        double[] polyF = getParamsCorrection(xF, yF, xGf, yGf, degX, degY, xp, yp);
        int sizex = (degX + 1) * (degX + 2) / 2;
        System.out.println("done.");
        return copyRef(polyF, 0, sizex-1);
    }

    public static double[] get_polynom2(double[] xF, double[] yF, double[] xGf, double[] yGf,
                                   int degX, int degY, double xp, double yp)
    {
        System.out.println("Obtaining correction polynomials... ");
        double[] polyF = getParamsCorrection(xF, yF, xGf, yGf, degX, degY, xp, yp);
        int sizex = (degX + 1) * (degX + 2) / 2;
        int sizey = (degY + 1) * (degY + 2) / 2;
        System.out.println("done.");
        return copyRef(polyF, sizex, sizex+sizey-1);
    }

    public static double[] minus(double[] mas1, double k){
        double[] mas = new double[mas1.length];
        for (int i=0; i< mas1.length; i++)
            mas[i]=mas1[i]-k;
        return mas;
    }

    public static double[] derive(double[] mas1, double k){
        double[] mas = new double[mas1.length];
        for (int i=0; i< mas1.length; i++)
            mas[i]=mas1[i]/k;
        return mas;
    }

    public static double mult(double[] v, double[] v1) {
        assert(v1.length == v.length);
        double res = 0;
        for(int i = v1.length-1; i >= 0; i--)
            res += v[i] * v1[i];
        return res;
    }

    public static double[] copyRef(double[] mas, int s1, int s2){
        assert(0 <= s1 && s1 <= s2 && s2 <= mas.length);
        double[] mm=new double[s2-s1+1];
        for (int i=0;i<mm.length;i++)
            mm[i]=mas[i+s1];
        return mm;
    }

    public static double[] getParamsCorrection(double[] x_corr, double[] y_corr, double[] x_dist, double[] y_dist, int degX, int degY, double xp, double yp)
    {
        int sizex = (degX+1)*(degX+2)/2;
        int sizey = (degY+1)*(degY+2)/2;
        int lenxy = x_corr.length;
        double[] x_dist_rad=minus(x_dist,xp);  double[] y_dist_rad = minus(y_dist,yp);
        double[] x_corr_rad=minus(x_corr,xp); double[] y_corr_rad = minus(y_corr,yp);
        double norm_xy = 0;
        for (int i = 0; i < lenxy; i++) norm_xy += x_dist_rad[i]*x_dist_rad[i] + y_dist_rad[i]*y_dist_rad[i];
        norm_xy = sqrt(norm_xy);
        x_dist_rad=derive(x_dist_rad, norm_xy); y_dist_rad  =derive(y_dist_rad, norm_xy);

        Matrix coefTermX=new Matrix(sizex, lenxy);
        Matrix coefTermY=new Matrix(sizey, lenxy);
        int idx = 0;
        for (int ii = degX; ii >= 0; ii--) {
            for (int j = 0; j <= ii; j++)  {
                for (int k = 0; k < lenxy; k++)
                    coefTermX.set(idx, k, pow(x_dist_rad[k], ii-j) * pow(y_dist_rad[k], j));
                idx++;
            }
        }
        idx = 0;
        for (int ii = degY; ii >= 0; ii--) {
            for (int j = 0; j <= ii; j++)  {
                for (int k = 0; k < lenxy; k++)
                    coefTermY.set(idx, k, pow(x_dist_rad[k], ii-j) * pow(y_dist_rad[k], j));
                idx++;
            }
        }

        int sizexy = sizex+sizey;
        Matrix coef_mat = new Matrix(sizexy, sizexy);      //rowref=row
        double[] m = new double[sizexy];
        for (int i=0;i<sizexy; i++)
            m[i]=0;
        int count_i = 0;
        for (int i = 0; i < sizex; i++) {
            int count_j = 0;
            for (int j = 0; j < sizex; j++) {
                coef_mat.set(count_i, count_j, mult(coefTermX.row(i), coefTermX.row(j)));
                count_j++; }
            m[count_i] = mult(x_corr_rad, coefTermX.row(i));
            count_i++; }

        count_i = 0;
        for (int i = 0; i < sizey; i++) {
            int count_j = 0;
            for (int j = 0; j < sizey; j++) {
                coef_mat.set(sizex+count_i, sizex+count_j, mult(coefTermY.row(i), coefTermY.row(j)));
                count_j++; }
            m[sizex+count_i] = mult(y_corr_rad, coefTermY.row(i));
            count_i++; 	}

        Matrix normalization_mat1 = new Matrix(sizexy, sizexy);
        Matrix inv_normalization_mat1 = new Matrix(sizexy, sizexy);
        for (int i = 0; i < sizexy; i++) {
            if (i < sizex)
                normalization_mat1.set(i,i, coef_mat.get(0, i));
            else
                normalization_mat1.set(i,i, coef_mat.get(sizex, i));
            inv_normalization_mat1.set(i,i, 1 / normalization_mat1.get(i,i));
        }

        Matrix normalized_coef_mat = coef_mat.mult(inv_normalization_mat1);
        Matrix normalization_mat2 = new Matrix(sizexy, sizexy);
        for (int i = 0; i < sizexy; i++) normalization_mat2.set(i, i,normalized_coef_mat.get(0, 0) / normalized_coef_mat.get(i, i));
        normalized_coef_mat =  normalization_mat2.mult(normalized_coef_mat);

        double[] paramsInv=new double[sizex+sizey];
        solveLU(normalized_coef_mat, normalization_mat2.mult_ret_vec(m), paramsInv);
        paramsInv = inv_normalization_mat1.mult_ret_vec(paramsInv);
        double[] denorm_paramsInv=new double[sizex+sizey];
        double[] denormX = copyRef(denorm_paramsInv,0, sizex-1);
        double[] denormY = copyRef(denorm_paramsInv, sizex, sizex+sizey-1);
        denormalization(denormX, denormY,  // copyRef is not const?!
                copyRef(paramsInv,0, sizex),
                copyRef(paramsInv, sizex, sizex+sizey-1),
                norm_xy, norm_xy, degX, degY);
        for (int i=0; i<sizex; i++)
            denorm_paramsInv[i]=denormX[i];
        for (int i=0; i<denormY.length;i++)
            denorm_paramsInv[i+sizex]=denormY[i];
        return denorm_paramsInv;
    }

    public static void denormalization(double[] denormX, double[] denormY, double[] normX, double[] normY,
                         double scaleX, double scaleY, int orderX, int orderY)
    {
        int idx = 0;
        for (int ii = orderX; ii >= 0; ii--) {
            for (int j = 1; j <= ii+1; j++) {
                denormX[idx] = normX[idx] / (pow(scaleX, ii-(j-1)) * pow(scaleY, j-1));
                idx++; }
        }
        idx = 0;
        for (int ii = orderY; ii >= 0; ii--) {
            for (int j = 1; j <= ii+1; j++) {
                denormY[idx] = normY[idx] / (pow(scaleX, ii-(j-1)) * pow(scaleY, j-1));
                idx++; }
        }
    }

    public static boolean solveLU(Matrix A, double[] B, double[] X)
    {
        //X = B;
        for (int i=0; i<X.length; i++)
            X[i]=B[i];
        return solveLU(A, X);
    }

    /// Replace X by A^{-1}X, by LU solver.
    public static boolean solveLU(Matrix A, double[] X)
    {
        assert(A.m_rows == A.m_cols);
        int	n = A.m_rows;
        double[] rowscale=new double[n]; // Implicit scaling of each row
        double[] permut =new double[n]; // Permutation of rows
        for (int i=0; i<permut.length; i++)
            permut[i]=0;

        // Get the implicit scaling information of each row
        for(int i=0; i< n; i++) {
            double max = 0.0;
            for (int j = 0; j < n; j++) {
                double tmp = abs(A.get(i, j));
                if (tmp > max)
                    max = tmp;
            }
            if (max == 0.0)
                return false;
            rowscale[i] = 1.0 / max;
        }
        // Perform the decomposition
        for(int k=0; k < n; k++) {
            // Search for largest pivot element
            double max = rowscale[k]*abs(A.get(k,k));
            int imax = k;
            for(int i=k+1; i < n; i++) {
                double tmp = rowscale[i]*abs(A.get(i,k));
                if(tmp > max) {
                    max = tmp;
                    imax = i;
                }
            }
            if(max == 0.0)
                return false;

            // Interchange rows if needed
            if(k != imax) {
                A.swapRows(k, imax);
                rowscale[imax] = rowscale[k]; // Scale of row k no longer needed
            }
            permut[k] = imax; // permut(k) was not initialized before
            double Akk = 1/A.get(k,k);
            for(int i=k+1; i < n; i++) {
                A.set(i,k, A.get(i,k)*Akk);
                double tmp = A.get(i,k); // Divide by pivot
                for (int j=k+1;j < n; j++) // Reduce the row
                    A.set(i,j,A.get(i,j)-tmp*A.get(k,j));
            }
        }
        // Forward substitution
        for (int k = 0; k < n; k++) {
            double sum = X[(int)permut[k]];
            X[(int)permut[k]] = X[k];
            for(int j = 0; j < k; j++)
                sum -= A.get(k,j)*X[j];
            X[k] = sum;
        }
        // Backward substitution
        for(int k = n - 1; k >= 0; k--) {
            double sum = X[k];
            for(int j=k+1; j < n; j++)
                sum -= A.get(k,j)*X[j];
            X[k] = sum/A.get(k,k);
        }
        return true;
    }

    public static void keypnts_circle(Image imgR, Image imgG, Image imgB, int scale, boolean clr) {
        int wiRB = imgR.xsize, heRB = imgR.ysize;
        int wiG = wiRB*scale, heG = heRB*scale;

        MyNum<Double> maxR, maxG, maxB;
        MyNum<Double> minR, minG, minB;
        maxR = new MyNum<Double>(0.0);
        maxG = new MyNum<Double>(0.0);
        maxB = new MyNum<Double>(0.0);
        minR = new MyNum<Double>(255.0);
        minG = new MyNum<Double>(255.0);
        minB = new MyNum<Double>(255.0);
        img_extremas(imgR, minR, maxR);
        img_extremas(imgG, minG, maxG);
        img_extremas(imgB, minB, maxB);
        double threR = 0.5*(maxR.get()-minR.get());
        double threG = 0.54*(maxG.get()-minG.get());
        double threB = 0.4*(maxB.get()-minB.get());

        Image imgbiR = new Image(wiRB, heRB, 255);
        Image imgbiG = new Image(wiG, heG, 255);
        Image imgbiB = new Image(wiRB, heRB, 255);

        binarization(imgbiR, imgbiG, imgbiB, imgR, imgG, imgB, threR, threG, threB);

        System.out.println("finding connected components... \n");
        ArrayList<CCStats> ccstatsR, ccstatsG, ccstatsB;
        ccstatsB= new ArrayList<>();
        ccstatsR= new ArrayList<>();
        ccstatsG= new ArrayList<>();
        CCStats.CC(ccstatsR, imgbiR, 'R');
        System.out.print(" number = " + ccstatsR.size() + " ");
        CCStats.CC(ccstatsG, imgbiG, 'G');
        System.out.print(" number = " + ccstatsG.size() + " ");
        CCStats.CC(ccstatsB, imgbiB, 'B');
        System.out.print(" number = " + ccstatsB.size() + " ");

        System.out.println("\nnumber of connected components per channel R, G, B = " + ccstatsR.size() + ", " + ccstatsG.size() +", " + ccstatsB.size());
        assert(ccstatsR.size() == ccstatsG.size() && ccstatsG.size() == ccstatsB.size());    //если неверно, программа завершится
        System.out.println("centers initialization for channels is done");

        System.out.println("\nMatching the centers... ");
        int ntaches = ccstatsG.size();
        xR = new double[ntaches]; yR = new double[ntaches]; rR = new double[ntaches];
        ones(xR); ones(yR); ones(rR);
        xB = new double[ntaches]; yB = new double[ntaches]; rB = new double[ntaches];
        ones(xB); ones(yB); ones(rB);
        xGr = new double[ntaches]; yGr = new double[ntaches]; rG = new double[ntaches];
        ones(xGr); ones(yGr); ones(rG);
        xGb = new double[ntaches]; yGb = new double[ntaches];
        ones(xGb); ones(yGb);
        for (int i = 0; i < ntaches; i++) {
            double xg = ccstatsG.get(i).centerX;
            double yg = ccstatsG.get(i).centerY;
            int idxR = CCStats.findMatch(xg, yg, ccstatsR, scale);
            int idxB = CCStats.findMatch(xg, yg, ccstatsB, scale);
            xGr[i] = xg; yGr[i] = yg;
            xGb[i] = xg; yGb[i] = yg;
            rG[i] = 0.5*(ccstatsG.get(i).radius1+ccstatsG.get(i).radius2);
            if (idxR != -1) {
                xR[i] = ccstatsR.get(idxR).centerX;
                yR[i] = ccstatsR.get(idxR).centerY;
                rR[i] = 0.5*(ccstatsR.get(idxR).radius1+ccstatsR.get(idxR).radius2); }
            if (idxB != -1) {
                xB[i] = ccstatsB.get(idxB).centerX;
                yB[i] = ccstatsB.get(idxB).centerY;
                rB[i] = 0.5*(ccstatsB.get(idxB).radius1+ccstatsB.get(idxB).radius2); }
        }
        System.out.println("done.");
        circle_redefine(imgR, imgG, imgB, scale, clr, true);
    }

    public static void circle_redefine(Image imgR, Image imgG, Image imgB, int scale, boolean clr, boolean green)
    {
        System.out.println("\nLMA center redefinition for the channels... \n");
        int ntaches = xR.length;
        for (int i = 0; i < ntaches; i++) {
            MyNum<Integer> x0R = new MyNum<>(0);
            MyNum<Integer> y0R = new MyNum<>(0);
            MyNum<Integer> x0G = new MyNum<>(0);
            MyNum<Integer> y0G = new MyNum<>(0);
            MyNum<Integer> x0B = new MyNum<>(0);
            MyNum<Integer> y0B = new MyNum<>(0);
            Image sub_imgR = takeSubImg(imgR, xR[i], yR[i], rR[i], x0R, y0R);   //левый нижний край
            Image sub_imgB = takeSubImg(imgB, xB[i], yB[i], rB[i], x0B, y0B);

            MyNum<Double> cxR, cyR, cxB, cyB;
            cxR = new MyNum<Double>(0.0);
            cyR = new MyNum<Double>(0.0);
            cxB = new MyNum<Double>(0.0);
            cyB = new MyNum<Double>(0.0);
            centerLMA (sub_imgR, clr, cxR, cyR);
            centerLMA (sub_imgB, clr, cxB, cyB);

            xR[i] = scale * (x0R.get() + cxR.get());
            yR[i] = scale * (y0R.get() + cyR.get());
            xB[i] = scale * (x0B.get() + cxB.get());
            yB[i] = scale * (y0B.get() + cyB.get());

            if (green) {
                Image sub_imgG = takeSubImg(imgG, xGr[i], yGr[i], rG[i], x0G, y0G);
                MyNum<Double> cxG, cyG;
                cxG = new MyNum<Double>(0.0);
                cyG = new MyNum<Double>(0.0);
                centerLMA (sub_imgG, clr, cxG, cyG);
                xGr[i] = x0G.get() + cxG.get();
                yGr[i] = y0G.get() + cyG.get();
                xGb[i] = x0G.get() + cxG.get();
                yGb[i] = y0G.get() + cyG.get();
            }
        }
    }

    public static double centerLMA(Image sub_img, boolean clr, MyNum<Double> centerX, MyNum<Double> centerY)
    {
        Image img_avg = average_image(sub_img);
        int w = sub_img.xsize;
        int h = sub_img.ysize;
        double cx = w/2, cy = h/2;
        MyNum<Double> radi = new MyNum<>(0.4*w);
        double[] P=new double[11];
        initial_tache(sub_img, P, radi, clr, cx, cy);
        //double[] trgData = trgtDataCalc(img_avg, P[3], P[4], radi.get()*2);
        //LMTacheC ellipseLMA = new LMTacheC(img_avg, P[3], P[4], radi.get()*2, clr, w, h);
        centerX.set(P[3]);
        centerY.set(P[4]);
        return 0;
    }

    public static double[] trgtDataCalc(Image img_avg, double cx, double cy, double delta) {
        int xbegin = (int)(cx+0.5-delta);
        int xend = (int)(cx+0.5+delta);
        int ybegin = (int)(cy+0.5-delta);
        int yend = (int)(cy+0.5+delta);
        int nerr = (yend-ybegin+1)*(xend-xbegin+1);
        double[] trgData = new double[nerr];
        for (int i = 0; i<trgData.length; i++)
            trgData[i]=0;
        int wi = img_avg.xsize;
        int he = img_avg.ysize;
        int idx = 0;
        for (int v = ybegin; v <= yend; v++) {
            for (int u = xbegin; u <= xend; u++) {
                if (u >= 1 && u <= wi-2 && v >= 1 && v <= he-2)
                    trgData[idx] = img_avg.data[u+v*img_avg.xsize];
                else
                    trgData[idx] = 255;
                idx++;
            }
        }
        return trgData;
    }

    public static int initial_tache(Image I, double[] h, MyNum<Double> rayon, boolean color, double x, double y) {
        int COL_IMA = I.xsize;
        int LIG_IMA = I.ysize;
        int j = (int)x;
        int i = (int)y;
        int d = (int)(2*rayon.get());
        if (2*d+1 > LIG_IMA)
            d=(LIG_IMA-1)/2;
        if(2*d+1 > COL_IMA)
            d=(COL_IMA-1)/2;
        if(i<d)
            i=d+1;
        if(i>LIG_IMA-1-d)
            i=LIG_IMA-2-d;
        if(j<d)
            j=d+1;
        if(j>COL_IMA-1-d)
            j=COL_IMA-2-d;
        int val_haut=0;
        int val_bas=255;
        for (int k = -d; k <= d; k++) {
            for (int l = -d;  l <= d; l++) {
                double lum = I.data[j+l+(i+k)*COL_IMA];
                if (lum > val_haut)
                    val_haut = (int)lum;
                else if (lum<val_bas)
                    val_bas=(int)lum;
            } }
        double seuil = 0;
        if (!color)
            seuil = val_bas + (val_haut - val_bas)/3 * 2;
        else
            seuil = val_bas + (val_haut - val_bas)/3;

        Matrix tab = new Matrix(2*d+1, 2*d+1);
        int label = 1;

        for (int k = -d+1; k <= d; k++){
            for(int l = -d+1; l <= d-1; l++){
                double lum= I.data[j+l+(i+k)*COL_IMA];
                if (lum < seuil) {
                    int imin=l;
                    int imax=l+1;
                    while( (I.data[j+imax+(i+k)*COL_IMA] <= seuil) && (imax <= d-1) ){
                        imax++; }
                    int vallab=0;
                    for(int m = imin; m <= imax; m++){
                        if (tab.get(k+d-1,m+d) != 0){
                            vallab = (int)tab.get(k+d-1,m+d);
                        }
                    }
                    if (vallab == 0){
                        vallab=label;
                        label++;
                    }
                    for(int m = imin; m <= imax; m++){
                        tab.set(k+d,m+d,vallab);
                    }
                    l=imax;
                }
            }
        }
        Matrix bary = new Matrix(label, 4);
        for(int k = -d; k <= d; k++){
            for(int l = -d; l <= d; l++){
                if(tab.get(k+d,l+d)!=0){
                    double lum= I.data[j+l+(i+k)*COL_IMA];
                    bary.set((int)tab.get(k+d,l+d),0,bary.get((int)tab.get(k+d,l+d),0)+(255-lum)*(j+l));
                    bary.set((int)tab.get(k+d,l+d),1,bary.get((int)tab.get(k+d,l+d),1)+(255-lum)*(i+k));
                    bary.set((int)tab.get(k+d,l+d),2, bary.get((int)tab.get(k+d,l+d),2)+(255-lum));
                    bary.set((int)tab.get(k+d,l+d),3, bary.get((int)tab.get(k+d,l+d),3)+1);
                }
            }
        }
        int distmin=100;
        int labelmin=0;
        for(int k = 1; k < label; k++){
            double dist = Math.sqrt( (bary.get(k,0)/bary.get(k,2)-x) * (bary.get(k,0)/bary.get(k,2)-x)+
                    (bary.get(k,1)/bary.get(k,2)-y) * (bary.get(k,1)/bary.get(k,2)-y));
            if(dist < distmin && bary.get(k,3) > 25 ){ /* 25 =  surface min*/
                distmin=(int)dist;
                labelmin=k;
            }
        }
        if(labelmin == 0) {
            System.out.println("not enough space (<=25 pixels)\n");   //слишком мало места
            return 1;
        }
        x=bary.get(labelmin,0)/bary.get(labelmin,2);
        y=bary.get(labelmin,1)/bary.get(labelmin,2);
        double sx2 = 0, sy2 = 0, sxy = 0, ss = 0;
        for(int k = -d; k <= d; k++){
            for(int l = -d; l <= d; l++){
                if(tab.get(k+d,l+d) == labelmin){
                    sx2+=(j+l-x)*(j+l-x);
                    sy2+=(i+k-y)*(i+k-y);
                    sxy+=(i+k-y)*(j+l-x);
                    ss++;
                }
            }
        }
        double lambda1 =  ((sx2+sy2)/ss + Math.sqrt(( (sx2+sy2)*(sx2+sy2)+4*(sxy*sxy-sx2*sy2)))/ss)/2.0;
        double lambda2 =  ((sx2*sy2-sxy*sxy)/(ss*ss))/lambda1;
        rayon.set(Math.sqrt(lambda1)*2);
        h[0] = 1.0/rayon.get(); 	                        /* lambda1 		*/
        h[1] = (Math.sqrt(lambda1/lambda2))/rayon.get(); 	/* lambda2		*/
        h[2] = Math.atan2(sx2/ss-lambda1, -sxy/ss);    	/* alpha  	 	*/
        h[3] = x;        	/* tu     		*/
        h[4] = y;        	/* tv      		*/
        h[5] = 0.25;       	/* радиус1  	*/
        h[6] = -2.0;      	/* откос     	*/
        h[7] = 0.25;       	/* радиус2   	*/
        h[8] = val_haut; 	/* val_haut   		*/
        h[9] = val_bas; 	/* val_bas   		*/
        h[10] = 1.0; 	/* position step  	*/
        return 0;
    }

    public static Image average_image(Image img) {
        int w = img.xsize;
        int h = img.ysize;
        Image img_avg = img.copy();
        for (int v = 1; v < h-1; v++) {
            for (int u = 1; u < w-1; u++) {
                double pix = img.data[u-1+(v-1)*w] + img.data[u+(v-1)*w] + img.data[u+1+(v-1)*w] +  //берём 9 пикселей, данный и 8 вокруг него
                        img.data[u-1+v*w] + img.data[u+v*w] + img.data[u+1+v*w] +
                                img.data[u-1+(v+1)*w] + img.data[u+(v+1)*w] + img.data[u+1+(v+1)*w] ;
                img_avg.data[u+v*w] = pix/9;
            }
        }
        return img_avg;
    }

    public static Image takeSubImg(Image IMG, double cx, double cy, double radi, MyNum<Integer> x0, MyNum<Integer> y0) //вырезает квадратную область с данным диском
    {
        int size = (int)(2.5 * radi);
        int x1 = (int)(cx - 0.5*size), x2 = (int)(cx + 0.5*size);  //в эту квадратную область залезет диск
        int y1 = (int)(cy - 0.5*size), y2 = (int)(cy + 0.5*size);
        x0.set(x1);
        y0.set(y1);
        if (y2-y1 != x2-x1)
            y2 = y1+x2-x1;
        Image img = new Image(x2-x1, y2-y1, 0);
        for (int i = 0; i < img.xsize; i++) {
            for (int j = 0; j < img.ysize; j++)
                if (x1+i >= 0 && y1+j >= 0 && x1+i<IMG.xsize && y1+j<IMG.ysize)
                    img.data[i+j*img.xsize] = IMG.data[x1+i+(y1+j)*IMG.xsize];
                else
                    img.data[i+j*img.xsize] = 255;
        }
        return img;
    }

    public static void ones(double[] mas){
        for (int i=0; i<mas.length; i++)
            mas[i]=1;
    }


    public static void binarization(Image imgbiR, Image imgbiG, Image imgbiB,
                                    Image imgR, Image imgG, Image imgB,
                      double threR, double threG, double threB)     //все тёмные участки закрашиваются чёрным
    {
        int wiRB = imgR.xsize;
        int heRB = imgR.ysize;
        int wiG = imgG.xsize;
        int heG = imgG.ysize;
        double red, blue, green;
        for (int i = 0; i < wiRB; i++) {
            for (int j = 0; j < heRB; j++) {
                red = imgR.data[i + j * wiRB];
                if (red <= threR) {
                    imgbiR.data[i + j * wiRB] = 0;        //lower than threshold, right
                }

                blue = imgB.data[i + j * wiRB];
                if (blue <= threB) imgbiB.data[i + j * wiRB] = 0;

                if (wiG == wiRB && heG == heRB) {
                    green = imgG.data[i+j*wiRB];
                    if (green <= threG) imgbiG.data[i+j*wiRB] = 0;	}
                else {
                    green = imgG.data[i*2+1+j*2*imgbiG.xsize];
                    if (green <= threG) imgbiG.data[i*2+1+j*2*imgbiG.xsize] = 0;
                    green = imgG.data[i*2+j*2*imgbiG.xsize];
                    if (green <= threG) imgbiG.data[i*2+j*2*imgbiG.xsize] = 0;
                    green = imgG.data[i*2+(j*2+1)*imgbiG.xsize];
                    if (green <= threG) imgbiG.data[i*2+(j*2+1)*imgbiG.xsize] = 0;
                    green = imgG.data[i*2+1+(j*2+1)*imgbiG.xsize];
                    if (green <= threG) imgbiG.data[i*2+1+(j*2+1)*imgbiG.xsize] = 0; }
            }
        }
    }

    public static void img_extremas(Image img, MyNum<Double> min, MyNum<Double> max) {
        min.set(255.0);
        max.set(0.0);
        for (int i = 4; i < img.xsize-4; i++) {
            for (int j = 4; j < img.ysize-4; j++) {
                double clr = img.data[j*img.xsize + i];
                if (clr > max.get()) max.set(clr);
                if (clr < min.get()) min.set(clr);
            }
        }
    }

    public static void raw2rgb(Image img_bayer, Image imgR, Image imgG, Image imgB){
        System.out.println("raw data processing...");
        int wiRB = imgR.xsize;
        int heRB = imgR.ysize;
        double red, blue, green;
        for (int i = 1; i < wiRB-1; i++) {
            for (int j = 1; j < heRB-1; j++) {
                red = img_bayer.data[i*2+j*2*img_bayer.xsize];
                imgR.data[i+j*wiRB] = red;

                blue = img_bayer.data[i*2+1+(j*2+1)*img_bayer.xsize];
                imgB.data[i+j*wiRB] = blue;

                green = img_bayer.data[i*2+1+j*2*img_bayer.xsize];
                imgG.data[i*2+1+j*2*imgG.xsize] = green;
                green = img_bayer.data[i*2+(j*2+1)*img_bayer.xsize];
                imgG.data[i*2+(j*2+1)*imgG.xsize] = green;
                imgG.data[i*2+j*2*imgG.xsize] = 0.25* (img_bayer.data[(i*2+1)+j*2*img_bayer.xsize] + img_bayer.data[i*2+(j*2+1)*img_bayer.xsize] +
                        img_bayer.data[(i*2-1)+j*2*img_bayer.xsize] + img_bayer.data[i*2+(j*2-1)*img_bayer.xsize]);
                imgG.data[i*2+1+(j*2+1)*imgG.xsize] = 0.25* (img_bayer.data[(i*2+1)+j*2*img_bayer.xsize] + img_bayer.data[i*2+(j*2+1)*img_bayer.xsize] +
                        img_bayer.data[(i*2+2)+(j*2+1)*img_bayer.xsize] + img_bayer.data[(i*2+1)+(j*2+2)*img_bayer.xsize]);
            }
        }
        System.out.println("done\n");
    }
}



