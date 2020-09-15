public class Matrix {
    public int m_rows; ///< Number of rows.
    public int m_cols; ///< Number of columns.
    public double[] p;

    public Matrix(int n, int m) {
        m_rows=n;
        m_cols=m;
        p = new double[n*m];
        for (int i=0; i<m*n; i++)
            p[i]=0;
    }

    public double get(int i, int j) {
        assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
        return p[i*m_cols + j];
    }

    public void set(int i, int j, double k){
        this.p[i*m_cols + j]= k;
    }

    public double[] col(int j){
        assert(j >= 0 && j < m_cols);
        double[] c=new double[m_rows];
        for(int i = 0; i < m_rows; i++) {
            c[i] = this.p[i*m_cols+j];
        }
        return c;
    }

    public void swapRows(int s1, int s2){
        assert(0 <= s1 && s1 < m_rows && 0 <= s2 && s2 < m_rows);
        double[] row1 = this.row(s1);
        double[] row2 = this.row(s2);
        for(int i = 0; i < m_cols; i++) {
            p[s1*m_cols+i] = row1[i];
            p[s2*m_cols+i] = row2[i];
        }
    }

    public double[] row(int j){
        assert(j >= 0 && j < m_rows);
        double[] c=new double[m_cols];
        for(int i = 0; i < m_cols; i++) {
            c[i] = this.p[j*m_cols+i];
        }
        return c;
    }

    public Matrix mult(Matrix b) {
        int m1 = this.m_rows;
        int n1 = this.m_cols;
        int m2 = b.m_rows;
        int n2 = b.m_cols;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix c = new Matrix(m1, n2);
        for (int i = 0; i < m1; i++)
            for (int j = 0; j < n2; j++)
                for (int k = 0; k < n1; k++)
                    c.p[i*n2+j] += this.p[i*n1+ k] * b.p[k*n2+j];
        return c;
    }

    public Matrix mult(double[] b) {
        Matrix m = new Matrix(b.length, 1);
        for (int i=0; i<b.length;i++) {
            m.p[i]=b[i];
        }
        return this.mult(m);
    }

    public double[] mult_ret_vec(double[] b) {
        Matrix mat = this.mult(b);
        double[] mas = new double[mat.m_rows*mat.m_cols];
        for (int i=0 ;i<mat.m_rows*mat.m_cols; i++)
            mas[i]=mat.p[i];
        return mas;
    }


    /// Tranposed of matrix
    public Matrix t() {
        Matrix t=new Matrix(m_cols, m_rows);
        for(int i = 0; i < t.m_rows; i++) {
            for(int j = 0; j < t.m_cols; j++) {
                t.p[i*t.m_cols+j]=this.p[j*t.m_rows+i];
            }
        }
        return t;
    }

}
