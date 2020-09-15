import java.util.ArrayList;
import java.util.Stack;

public class CCStats {
    public int nPoints;
    public double centerX, centerY, radius1, radius2;
    public int perimeter; // need to calculate compactness measure of shape to eliminate noise

    public static int white_neighbors(Pixel p, Image img) {      //вокруг есть хотя бы один белый пиксель(возможный край диска)
        int nn = 0;
        double color;
        Pixel p1= new Pixel(p.x-1, p.y);
        if (p1.x >= 0 && p1.x < img.xsize && p1.y >= 0 && p1.y < img.ysize) {
            if (img.data[p1.x+p1.y*img.xsize] == 255)
                nn++;
        } else nn++;                       //край

        Pixel p2 = new Pixel(p.x+1, p.y);
        if (p2.x >= 0 && p2.x < img.xsize && p2.y >= 0 && p2.y < img.ysize) {
            if (img.data[p2.x+p2.y*img.xsize] == 255)
                nn++;
        } else nn++;

        Pixel p3 = new Pixel(p.x, p.y-1);
        if (p3.x >= 0 && p3.x < img.xsize && p3.y >= 0 && p3.y < img.ysize) {
            if (img.data[p3.x + p3.y * img.xsize] == 255)
                nn++;
        } else nn++;

        Pixel p4 = new Pixel(p.x, p.y+1);
        if (p4.x >= 0 && p4.x < img.xsize && p4.y >= 0 && p4.y < img.ysize) {
            if (img.data[p4.x+p4.y*img.xsize] == 255)
                nn++;
        } else nn++;

        if (nn > 0) nn = 1;
        return nn;
    }

    public static void extract_CCStats(ArrayList<Pixel> cc, CCStats stats, Image img) {
        stats.nPoints = cc.size();
        stats.perimeter = 0;
        double meanX = 0, meanY = 0;
        double minX = img.xsize, minY = img.ysize, maxX = 0, maxY = 0;
        for (int i = 0; i < cc.size(); i++) {
            meanX += cc.get(i).x;
            meanY += cc.get(i).y;
            if (cc.get(i).x > maxX) maxX = cc.get(i).x;
            if (cc.get(i).x < minX) minX = cc.get(i).x;
            if (cc.get(i).y > maxY) maxY = cc.get(i).y;
            if (cc.get(i).y < minY) minY = cc.get(i).y;
            Pixel tmp = new Pixel(cc.get(i).x, cc.get(i).y);
            stats.perimeter += white_neighbors(tmp, img);      //по факту считается длина окружности, то есть суммируются все граничные пиксели
        }
        stats.centerX = meanX / cc.size();
        stats.centerY = meanY / cc.size();
        stats.radius1 = 0.5 * (maxX - minX);
        stats.radius2 = 0.5 * (maxY - minY);
    }

    public static int extract_cc_(Pixel p, ArrayList<Pixel> cc, Image img) {  //в массив пикселей засовываются чёрные пиксели, то есть область диска
        Stack<Pixel> s= new Stack<>();
        double color; int px, py;
        if (p.x >= 0 && p.x < img.xsize && p.y >= 0 && p.y < img.ysize) {
            color = img.data[p.x + p.y * img.xsize];
            px=p.x;py=p.y;
        } else
            return 0;
        while (color == 0 || !s.empty()) {
            if (color == 0) {
                cc.add(p);
                color=255;
                img.data[p.x + p.y * img.xsize]=255;
                    s.push(new Pixel(p.x + 1, p.y));
                    s.push(new Pixel(p.x - 1, p.y));
                    s.push(new Pixel(p.x, p.y+1));
                    s.push(new Pixel(p.x, p.y-1));
                    px=p.x;py=p.y;
            }
            if (!s.empty()) {
                p = s.pop();
                if (p.x >= 0 && p.x < img.xsize && p.y >= 0 && p.y < img.ysize)
                    color = img.data[p.x + p.y * img.xsize];
                else {
                    color = 255;
                    img.data[px + py * img.xsize]=255;
                }
            }
        }
        return cc.size();
    }

    public static void CC(ArrayList<CCStats> ccstats, Image imgbi, char channel) {
        Image img_copy = imgbi.copy();
        System.out.println("\nchannel " + channel + " : ");
        double meansize = 0;
        for (int i = 0; i < imgbi.xsize; i++) {
            for (int j = 0; j < imgbi.ysize; j++) {
                ArrayList<Pixel> ccC = new ArrayList<>();
                CCStats stats = new CCStats();
                Pixel tmp = new Pixel(i, j);
                int npix = extract_cc_(tmp, ccC, imgbi);
                if (npix > 180) {            //если это диск
                    extract_CCStats(ccC, stats, img_copy);       //запишет в stats центр диска
                    double compactness = 4 * Math.PI * stats.nPoints / (stats.perimeter * stats.perimeter);     //коэффициент сжатия эллипса
                    if (Math.min(stats.radius1, stats.radius2) > 8 && compactness < 1.3 && compactness > 0.7) {
                        ccstats.add(stats);
                        meansize += stats.nPoints;
                    }
                }
            }
        }
    }

    // finds corresponding circle center match among two channels
    public static int findMatch(double xg, double yg, ArrayList<CCStats> trgstats, double scale) {
        int minidx = -1;
        double mindist = 5; // look for center in radius of 10 pixels
        for (int i = 0; i < trgstats.size(); i++) {
            double eucdist = Math.sqrt((scale * trgstats.get(i).centerX - xg) * (scale * trgstats.get(i).centerX - xg) + (scale * trgstats.get(i).centerY - yg) * (scale * trgstats.get(i).centerY - yg));
            if (eucdist <= mindist) {
                mindist = eucdist;
                minidx = i;
            }
        }
        return minidx;
    }
}

 class Pixel {
     public int x, y;
     Pixel(int X, int Y) {
         x = X;
         y = Y;
     }
 }