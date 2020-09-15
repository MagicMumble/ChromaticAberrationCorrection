import java.io.*;

import static sun.nio.ch.IOStatus.EOF;

public class Image {
    public double[] data;
    public int xsize, ysize;

    public Image(int x, int y, double fill){
        if (x == 0 || y == 0)
            throw new Error("invalid image size");
        data=new double[x*y];
        xsize = x;
        ysize = y;
        for (int i=0; i<x*y; i++)
            this.data[i]=fill;
    }

    public Image copy() {
        Image im = new Image(this.xsize, this.ysize, 0);
        if (xsize * ysize >= 0)
            System.arraycopy(this.data, 0, im.data, 0, xsize * ysize);
        return im;
    }

    public static Image read_pgm(String name) throws IOException {

        boolean bin = false;
        int xsize, ysize, depth, x,y;
        Image im;
        File file = new File(name);
        FileInputStream reader = new FileInputStream(file);
        int c;
        if (reader.read() != 'P')
            throw new Error("not a PGM file!");
        if ((c=reader.read()) == '2')
            bin = false;
        else if ( c == '5')
            bin = true;
        else
            throw new Error("not a PGM file!");
        skip_whites_and_comments(reader);
        xsize = get_num(reader, -1);            /* X size */
        int ch = skip_whites_and_comments(reader);
        ysize = get_num(reader, ch);            /* Y size */
        ch = skip_whites_and_comments(reader);
        depth = get_num(reader, ch);
        if(depth==0)
                throw new Error("Warning: depth=0, probably invalid PGM file\n");
            im = new Image(xsize, ysize, 0);
            for(y=0;y<ysize;y++)
                for(x=0;x<xsize;x++) {
                    im.data[x + y * xsize] = bin ? (double) reader.read()
                            : (double) get_num(reader, -1);
                }
        return im;
    }


    private static int get_num(FileInputStream f, int ch) throws IOException {
        int num, c;
        while((c=f.read())==' ');
        if(!Character.isDigit(c))
            throw new Error("corrupted PGM file.");
        if (ch!=-1) {
            num = (ch - '0');
            num = 10 * num + c - '0';
        } else {
            num = (c - '0');
        }
        while( Character.isDigit(c=f.read()))
            num = 10 * num + c - '0';
        return num;
    }

    private static int skip_whites_and_comments(FileInputStream f) throws IOException {
        int c;
        do
        {
            while((c=f.read())==' '); /* skip spaces */
            if(c=='#') { /* skip comments */
                while (c != '\n' && c != '\r' && c != EOF)
                    c = f.read();
                c = f.read();
            }
        }
        while( c == '#' || c==' ' );
        return c;
    }
}
