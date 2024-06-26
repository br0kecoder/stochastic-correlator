import java.util.*;
public class SCCCorrelatorFormattedText
{
    private static double SNToProb(int[] x){
        int n=x.length;
        int c=0;
        for (int v:x) c+=v;
        return c*1.0/n;
    }
    private static int[] probToSN(double prob, int n){
        int[] seed1={0,0,0,1};
        int[][] rng1=lfsr(seed1,2,3);
        int c=(int)Math.round(prob*n);
        int[] x_bin=decToBin(c,4);
        int[] x_stoch=comparator(rng1,x_bin);
        return x_stoch;
    }
    private static int[] decToBin(int dec, int n){
        int[] bin=new int[n];
        while (dec>0){
            bin[n-1]=dec%2;
            dec/=2;
            n--;
        }
        return bin;
    }
    private static int[][] lfsr(int[] seed, int tap1, int tap2){
        int[][] seq=new int[16][4];
        for (int i=0;i<16;i++){
            //System.out.println(i);
            for (int j=0;j<4;j++){
                seq[i][j]=seed[j];
                //System.out.print(seq[i][j]+" ");
            }
            //System.out.println();
            int temp=seed[0];
            seed[0]=seed[tap1]^seed[tap2];
            seed[3]=seed[2];
            seed[2]=seed[1];
            seed[1]=temp;
            
            //System.out.print(seed[tap1]+" "+seed[tap2]+" "+(seed[tap1]^seed[tap2]));
            //System.out.println();
        }
        return seq;
    }
    private static double scc(int[] x, int[] y, int n){
        int a=0, b=0, c=0, d=0;
        double ans=0;
        for (int i=0;i<n;i++){
            if (x[i]==1&&y[i]==1) a++;
            else if(x[i]==0&&y[i]==0) d++;
            else if (x[i]==1&&y[i]==0) b++;
            else c++;
        }
        System.out.println("\n1-1 overlap="+a+"\n1-0 overlap="+b+"\n0-1 overlap="+c+"\n0-0 overlap="+d);
        if (a*d>b*c){
            if (n*Math.min(a+b,a+c)-(a+b)*(a+c)==0){
                if ((a*d-b*c)*1.0==0) ans=0;
                else if ((a*d-b*c)*1.0<0) ans=-1;
                else ans=1;
            }
            
            else ans=((a*d-b*c)*1.0)/(n*Math.min(a+b,a+c)-(a+b)*(a+c));
        }
        else{
            if ((a+b)*(a+c)-n*Math.max(a-d,0)==0){
                if ((a*d-b*c)*1.0==0) ans=0;
                else if ((a*d-b*c)*1.0<0) ans=-1;
                else ans=1;
            }
            ans=((a*d-b*c)*1.0)/((a+b)*(a+c)-n*Math.max(a-d,0));
        }
        return ans;
    }
    private static int[] comparator(int[][] rng, int[] bits){
        int sn[]=new int[16];
        for (int i=0;i<16;i++){
            int[] cur=rng[i];
            int inbit=0, rngbit=0, deg=1;
            for (int j=3;j>=0;j--){
                rngbit+=cur[j]*deg;
                inbit+=bits[j]*deg;
                deg*=2;
            }
            //System.out.println(rngbit+" "+inbit);
            if (inbit>rngbit) sn[i]=1;
            else sn[i]=0;
        }
        return sn;
    }
    private static void displaySN(int[] sn){
        for (int x:sn) System.out.print(x+" ");
        System.out.println();
    }
    private static int[][] invert(int[][] rng){
        for (int i=0;i<rng.length;i++){
            for (int j=0;j<rng[0].length;j++){
                if (rng[i][j]==1) rng[i][j]=0;
                else rng[i][j]=1;
            }
        }
        return rng;
    }
    private static int[] zeroSCC(int[] x, int[] y, int n){
        int x1=0, y1=0;
        for (int z:x) x1+=z;
        for (int z:y) y1+=z;
        System.out.println("x1="+x1+" y1="+y1);
        int a=(int)Math.round((x1*y1*1.0)/n), b=x1-a, c=y1-a, d=n-a-b-c;
        System.out.println("a=(x1*y1)/n="+a+" b=x1-a="+b+" c=y1-a="+c+" d=n-a-b-c="+d);
        //Arrays.fill(y,0);
        for (int i=0;i<n;i++){
            if (x[i]==1){
                if (y[i]==1){
                    if (a>0) a--;
                    else{
                        y[i]=0;
                        b--;
                    }
                }
                else{
                    if (b>0) b--;
                    else{
                        y[i]=1;
                        a--;
                    }
                }
            }
            else{
                if (y[i]==1){
                    if (c>0) c--;
                    else{
                        y[i]=0;
                        d--;
                    }
                }
                else{
                    if (d>0) d--;
                    else{
                        y[i]=1;
                        c--;
                    }
                }
            }
        }
        return y;
        
    }
    private static int[] anySCC(int[] x, int[] y, int n, double s){
        int x1=0, y1=0;
        int y_new[]=new int[y.length];
        for (int z:x) x1+=z;
        for (int z:y) y1+=z;
        //System.out.println("x1="+x1+" y1="+y1);
        int a=0;
        if (s>=0){
            if (y1>x1) a=(int)Math.round((x1*(s*n-s*y1+y1)*1.0)/n);
            else a=(int)Math.round((y1*(s*n-s*y1+y1)*1.0)/n);
        }
        else{
            if (x1+y1<=n) a=(int)Math.round((x1*y1*(s+1)*1.0)/n);
            else a=(int)Math.round((x1*y1*(s+1)*1.0-n*s*(x1+y1-n)*1.0)/n);
        }
        //(int)Math.round((x1*y1*1.0)/n);
        int b=x1-a, c=y1-a, d=n-a-b-c;
        if (d<0){
            a++;
            b--;
            c--;
            d=n-a-b-c;
        }
        if (a<0){
            b+=Math.abs(a);
            a=0;
        }
        else if (b<0){
            a+=Math.abs(b);
            b=0;
        }
        if (c<0){
            d+=Math.abs(c);
            c=0;
        }
        System.out.println("Required 1-1="+a+" 1-0=x1-a="+b+" 0-1=y1-a="+c+" 0-0=n-a-b-c="+d);
        //Arrays.fill(y,0);
        for (int i=n-1;i>=0;i--){
            int temp=y[i];
            if (x[i]==1){
                if (y[i]==1){
                    System.out.print("1-1: a="+a+"=>");
                    if (a>0) {a--; System.out.print("a="+a);}
                    else{
                        y[i]=0;
                        b--;
                        System.out.print("b="+b);
                    }
                }
                else{
                    System.out.print("1-0: b="+b+"=>");
                    if (b>0) {b--; System.out.print("b="+b);}
                    else{
                        y[i]=1;
                        a--;
                        System.out.print("a="+a);
                    }
                }
            }
            else{
                
                if (y[i]==1){
                    System.out.print("0-1: c="+c+"=>");
                    if (c>0) {c--; System.out.print("c="+c);}
                    else{
                        y[i]=0;
                        d--;
                        System.out.print("d="+d);
                    }
                }
                else{
                    System.out.print("0-0: d="+d+"=>");
                    if (d>0) {d--; System.out.print("d="+d);}
                    else{
                        y[i]=1;
                        c--;
                        System.out.print("c="+c);
                    }
                }
            }
            y_new[i]=y[i];
            for (int space=1;space<=(i+1)*2;space++) System.out.print(" ");
            System.out.println(y_new[i]);
            //displaySN(y);
            y[i]=temp;
        }
        int y2=0;
        for (int z:y_new) y2+=z;
        //System.out.println("Final y value: "+y2);
        return y_new;
        
    }
    private static int[] synch(int[] x, int[] y, int dep){ //only changes value of y
        int[] y_new=new int[y.length];
        int save0=0, save1=0, yval=0, yval_new=0;
        for (int i=0;i<y.length;i++){
            if (x[i]==y[i]) y_new[i]=y[i];
            else if (y[i]==1){
                if (save0>0){
                    save0--;
                    y_new[i]=0;
                }
                else if (save1+save0<dep){
                    save1++;
                    y_new[i]=0;
                }
                else y_new[i]=1;
            }
            else{
                if (save1>0){
                    save1--;
                    y_new[i]=1;
                }
                else if (save0+save1<dep){
                    save0++;
                    y_new[i]=1;
                }
                else y_new[i]=0;
            }
            yval+=y[i];
            yval_new+=y_new[i];
        }
        double bias=(yval-yval_new)*1.0/y.length;
        System.out.println("Bias="+bias);
        return y_new;
    }
    private static int[] desynch(int[] x, int[] y, int dep){ //only changes value of y
        int[] y_new=new int[y.length];
        int save0=0, save1=0, yval=0, yval_new=0;
        for (int i=0;i<y.length;i++){
            if (x[i]!=y[i]) y_new[i]=y[i];
            else if (y[i]==1){
                if (save0>0){
                    save0--;
                    y_new[i]=0;
                }
                else if (save1+save0<dep){
                    save1++;
                    y_new[i]=0;
                }
                else y_new[i]=1;
            }
            else{
                if (save1>0){
                    save1--;
                    y_new[i]=1;
                }
                else if (save0+save1<dep){
                    save0++;
                    y_new[i]=1;
                }
                else y_new[i]=0;
            }
            yval+=y[i];
            yval_new+=y_new[i];
        }
        double bias=(yval-yval_new)*1.0/y.length;
        System.out.println("Bias="+bias);
        return y_new;
    }
    private static double getBias(int[] x, int[] y){
        int xval=0, yval=0;
        for (int i=0;i<x.length;i++){
            xval+=x[i];
            yval+=y[i];
        }
        double b=(xval-yval)*1.0/x.length;
        return b;
    }
    private static int[][] synch_xy(int[] x, int[] y, int dep){ //changes value of both x and y
        int state=1;
        int[] x_new=new int[x.length];
        int[] y_new=new int[y.length];
        int xy[][]=new int[2][x.length];
        for (int i=0;i<x.length;i++){
            if (state==1){
                if (x[i]==y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==1 && y[i]==0){
                    x_new[i]=0;
                    y_new[i]=0;
                    state=0;
                }
                else{
                    x_new[i]=0;
                    y_new[i]=0;
                    state=2;
                }
            }
            else if (state==0){
                if (x[i]==y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==1 && y[i]==0){
                    x_new[i]=1;
                    y_new[i]=0;
                    
                }
                else{
                    x_new[i]=1;
                    y_new[i]=1;
                    state=1;
                }
            }
            else{
                if (x[i]==y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==0 && y[i]==1){
                    x_new[i]=0;
                    y_new[i]=1;
                }
                else{
                    x_new[i]=1;
                    y_new[i]=1;
                    state=1;
                }
            }
            xy[0][i]=x_new[i];
            xy[1][i]=y_new[i];
        }
        return xy;
    }
    private static int[][] desynch_xy(int[] x, int[] y, int dep){
        int state=0, n=x.length;
        int x_new[]=new int[n];
        int y_new[]=new int[n];
        int xy[][]=new int[2][n];
        for (int i=0;i<n;i++){
            if (state==0){
                if (x[i]!=y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==0 && y[i]==0){
                    x_new[i]=0;
                    y_new[i]=0;
                }
                else{
                    x_new[i]=0;
                    y_new[i]=1;
                    state=1;
                }
            }
            else if (state==1){
                if (x[i]!=y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==1&&y[i]==1){
                    x_new[i]=1;
                    y_new[i]=1;
                }
                else{
                    x_new[i]=1;
                    y_new[i]=0;
                    state=2;
                }
            }
            else if (state==2){
                if (x[i]!=y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==0&&y[i]==0){
                    x_new[i]=0;
                    y_new[i]=0;
                }
                else{
                    x_new[i]=1;
                    y_new[i]=0;
                    state=3;
                }
            }
            else{
                if (x[i]!=y[i]){
                    x_new[i]=x[i];
                    y_new[i]=y[i];
                }
                else if (x[i]==1&&y[i]==1){
                    x_new[i]=1;
                    y_new[i]=1;
                }
                else{
                    x_new[i]=0;
                    y_new[i]=1;
                    state=0;
                }
            }
            xy[0][i]=x_new[i];
            xy[1][i]=y_new[i];
        }
        return xy;
    }
    private static int[][] decorrelate(int[] x, int[] y, int dep){
        int n=x.length;
        int[][] xy=new int[2][n];
        
        int[] seed1={0,0,0,1};
        int[][] rng1=lfsr(seed1,2,3);        
        int[] seed2={0,0,0,1};
        int[][] rng2=lfsr(seed2,0,3);
        
        int[] sb1_bin=decToBin(8,4);
        int[] sb2_bin=decToBin(8,4);
        int[] sbx=comparator(rng1,sb1_bin);
        int[] sby=comparator(rng2,sb2_bin);
        
        int x_new[]=new int[n];
        int y_new[]=new int[n];
        
        int x0=0, x1=0, y0=0, y1=0;
        for (int i=0;i<n;i++){
            //System.out.print(x[i]);
            if (x[i]==sbx[i]&&y[i]==sby[i]){
               
                xy[0][i]=x[i];
                xy[1][i]=y[i];
                continue;
            }
            
            if (x[i]==0&&sbx[i]==1){
                if (x1>0){
                    //System.out.println("Y");
                    x_new[i]=1;
                    x1--;
                    
                }
                else if (x0+x1<dep){
                    //System.out.println("Y"+x[i]+" "+i);
                    x_new[i]=1;
                    x0++;
                    //System.out.println(x[i]);
                }
                else x[i]=0;
            }
            else if (x[i]==1&&sbx[i]==0){
                if (x0>0){
                    //System.out.println("Y");
                    x_new[i]=0;
                    x0--;
                }
                else if (x0+x1<dep){
                    //System.out.println("Y");
                    x_new[i]=0;
                    x1++;
                }
                else x_new[i]=1;
            }
            xy[0][i]=x_new[i];
            
            if (y[i]==0&&sby[i]==1){
                if (y1>0){
                    y_new[i]=1;
                    y1--;
                }
                else if (y0+y1<dep){
                    y_new[i]=1;
                    y0++;
                }
                else y[i]=0;
            }
            else if (y[i]==1&&sby[i]==0){
                if (y0>0){
                    y_new[i]=0;
                    y0--;
                }
                else if (y0+y1<dep){
                    y_new[i]=0;
                    y1++;
                }
                else y_new[i]=1;
            }
            xy[1][i]=y_new[i];
        }
        return xy;
    }
    public static void synchDesynch(){
        int[] seed1={0,0,0,1};
        int[][] rng1=lfsr(seed1,2,3);        
        int[] seed2={0,0,0,1};
        int[][] rng2=lfsr(seed2,0,3);
        
        int count=0;
        double totXBias=0, totYBias=0, totSCC=0, inSCC=0, totXAbsBias=0, totYAbsBias=0, mse_x=0, mse_y=0, abs_scc_dev=0, max_x_bias=0, max_y_bias=0, max_scc_bias=0;
        for (int x1=2;x1<=15;x1++){
            for (int y1=2;y1<=15;y1++){
                count++;
                
                int[] x_bin=decToBin(x1,4);
                int[] y_bin=decToBin(y1,4);
                int[] x_stoch=comparator(rng1,x_bin);
                int[] y_stoch=comparator(rng2,y_bin); 
                
                
                /*int[][] xy_synch=synch_xy(x_stoch,y_stoch,1); //synchronizer
                int[] x_synch=xy_synch[0];
                int[] y_synch=xy_synch[1];
                totXBias+=getBias(x_stoch,x_synch);
                totYBias+=getBias(y_stoch,y_synch);
                totSCC+=scc(x_synch,y_synch,16);
                inSCC+=scc(x_stoch,y_stoch,16);
                totXAbsBias+=Math.abs(getBias(x_stoch,x_synch));
                totYAbsBias+=Math.abs(getBias(y_stoch,y_synch));
                mse_x+=getBias(x_stoch,x_synch)*getBias(x_stoch,x_synch);
                mse_y+=getBias(y_stoch,y_synch)*getBias(y_stoch,y_synch);
                abs_scc_dev+=Math.abs(1-scc(x_synch,y_synch,16));
                if (Math.abs(getBias(x_stoch,x_synch))>Math.abs(max_x_bias))
                    max_x_bias=getBias(x_stoch,x_synch);
                if (Math.abs(getBias(y_stoch,y_synch))>Math.abs(max_y_bias))
                    max_y_bias=getBias(y_stoch,y_synch);
                if (Math.abs(scc(x_synch,y_synch,16)-1)>Math.abs(max_scc_bias))
                    max_scc_bias=scc(x_synch,y_synch,16)-1;*/
                
                //System.out.println(scc(x_synch,y_synch,16));
                
                
                /*int[][] xy_desynch=desynch_xy(x_stoch,y_stoch,1); //desynchronizer
                int[] x_desynch=xy_desynch[0];
                int[] y_desynch=xy_desynch[1];
                totXBias+=getBias(x_stoch,x_desynch);
                totYBias+=getBias(y_stoch,y_desynch);
                totSCC+=scc(x_desynch,y_desynch,16);
                inSCC+=scc(x_stoch,y_stoch,16);
                totXAbsBias+=Math.abs(getBias(x_stoch,x_desynch));
                totYAbsBias+=Math.abs(getBias(y_stoch,y_desynch));
                mse_x+=getBias(x_stoch,x_desynch)*getBias(x_stoch,x_desynch);
                mse_y+=getBias(y_stoch,y_desynch)*getBias(y_stoch,y_desynch);
                abs_scc_dev+=Math.abs(-1-scc(x_desynch,y_desynch,16));
                if (Math.abs(getBias(x_stoch,x_desynch))>Math.abs(max_x_bias))
                    max_x_bias=getBias(x_stoch,x_desynch);
                if (Math.abs(getBias(y_stoch,y_desynch))>Math.abs(max_y_bias))
                    max_y_bias=getBias(y_stoch,y_desynch);
                if (Math.abs(scc(x_desynch,y_desynch,16)+1)>Math.abs(max_scc_bias))
                    max_scc_bias=scc(x_desynch,y_desynch,16)+1;*/
                
                int[][] xy_dec=decorrelate(x_stoch,y_stoch,4); //decorrelator
                int[] x_dec=xy_dec[0];
                int[] y_dec=xy_dec[1];
                totXBias+=getBias(x_stoch,x_dec);
                totYBias+=getBias(y_stoch,y_dec);
                totSCC+=scc(x_dec,y_dec,16);
                inSCC+=scc(x_stoch,y_stoch,16);
                totXAbsBias+=Math.abs(getBias(x_stoch,x_dec));
                totYAbsBias+=Math.abs(getBias(y_stoch,y_dec));
                mse_x+=getBias(x_stoch,x_dec)*getBias(x_stoch,x_dec);
                mse_y+=getBias(y_stoch,y_dec)*getBias(y_stoch,y_dec);
                abs_scc_dev+=Math.abs(0-scc(x_dec,y_dec,16));
                if (Math.abs(getBias(x_stoch,x_dec))>Math.abs(max_x_bias))
                    max_x_bias=getBias(x_stoch,x_dec);
                if (Math.abs(getBias(y_stoch,y_dec))>Math.abs(max_y_bias))
                    max_y_bias=getBias(y_stoch,y_dec);
                if (Math.abs(0-scc(x_dec,y_dec,16))>Math.abs(max_scc_bias))
                    max_scc_bias=0-scc(x_dec,y_dec,16);
                    
                /*double s=-0.8;
                
                int[] y_anySCC=anySCC(x_stoch,y_stoch,16,s);
                totYBias+=getBias(y_stoch,y_anySCC);
                totSCC+=scc(x_stoch,y_anySCC,16);
                inSCC+=scc(x_stoch,y_stoch,16);
                totYAbsBias+=Math.abs(getBias(y_stoch,y_anySCC));
                mse_y+=getBias(y_stoch,y_anySCC)*getBias(y_stoch,y_anySCC);*/
            }
        }
        System.out.println("SCC=0.0, Decorrelator:");
        System.out.println("Average X Bias: "+totXBias/count);
        System.out.println("Average Y Bias: "+totYBias/count);
        System.out.println("Average absolute X Bias: "+totXAbsBias/count);
        System.out.println("Average absolute Y Bias: "+totYAbsBias/count);
        System.out.println("X MSE: "+mse_x/count);
        System.out.println("Y MSE: "+mse_y/count);
        System.out.println("Average output SCC: "+totSCC/count);
        System.out.println("Average input SCC: "+inSCC/count);
        System.out.println("Average absolute SCC deviation: "+abs_scc_dev/count);
        System.out.println("Max X Bias: "+max_x_bias);
        System.out.println("Max Y Bias: "+max_y_bias);
        System.out.println("Max SCC Bias: "+max_scc_bias);
        /*int[] xin={1,0,0,1,1,1,1,0,0,1,0,0,0,0,1,0};
        int[] yin={0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,1};
        int[][] xyin=decorrelate(xin,yin,1);
        displaySN(xin);
        displaySN(yin);
        displaySN(xyin[0]);
        displaySN(xyin[1]);
        System.out.println(scc(xyin[0],xyin[1],16));        */
    }
    private static int[] multiply(int[] x, int[] y){ //AND gate
        int[] y_uncor=anySCC(x,y,16,0);
        int[] xy_and=new int[x.length];
        for (int i=0;i<x.length;i++){
            xy_and[i]=x[i]&y_uncor[i];
        }
        return xy_and;
    }
    private static int[] absDiff(int[] x, int[] y){ //XOR gate
        int[] y_posSC=anySCC(x,y,16,1);
        int[] xy_xor=new int[x.length];
        for (int i=0;i<x.length;i++){
            xy_xor[i]=x[i]^y_posSC[i];
        }
        return xy_xor;
    }
    public static void mlc(){
        Scanner sc=new Scanner(System.in);
        double px[]=new double[5];
        for (int i=1;i<5;i++){
            System.out.print("Enter px"+(i)+": ");
            px[i]=sc.nextDouble();
        }
        
        double true_val=px[4]*Math.abs(px[1]*px[2]-px[2]*px[3]);
        System.out.println("Actual value of ckt: "+true_val);
        
        int[] x1=probToSN(px[1], 16);
        int[] x2=probToSN(px[2], 16);
        int[] x3=probToSN(px[3], 16);
        int[] x4=probToSN(px[4], 16);
        
        int[] y1=multiply(x1,x2);
        int[] y2=multiply(x2,x3);
        
        int[] z1=absDiff(y1,y2);
        int[] z2=multiply(z1,x4);
        
        System.out.println("Value after operation: "+SNToProb(z2));
    }
    public static void mlc2(){
        Scanner sc=new Scanner(System.in);
        double px[]=new double[3];
        for (int i=1;i<3;i++){
            System.out.print("Enter px"+(i)+": ");
            px[i]=sc.nextDouble();
        }
        
        double true_val=px[1]*px[1]*px[1]*px[2]*px[2]*px[2];
        System.out.println("Actual value of ckt: "+true_val);
        
        int[] x1=probToSN(px[1], 16);
        int[] x2=probToSN(px[2], 16);
        int[] x3=probToSN(px[3], 16);
        int[] x4=probToSN(px[4], 16);
        
        int[] y1=multiply(x1,x2);
        int[] y2=multiply(x2,x3);
        
        int[] z1=absDiff(y1,y2);
        int[] z2=multiply(z1,x4);
        
        System.out.println("Value after operation: "+SNToProb(z2));
    }
    public static void main_2(String args[]){
        Scanner sc=new Scanner(System.in);
        int[] seed1={0,0,0,1};
        int[][] rng1=lfsr(seed1,2,3);        
        int[] seed2={0,0,0,1};
        int[][] rng2=lfsr(seed2,0,3);
        
        while (true){
        System.out.println("*********************************");
        System.out.print("Enter x in decimal: ");
        int x1=sc.nextInt();
        System.out.print("Enter y in decimal: ");
        int y1=sc.nextInt();
        System.out.print("Enter save depth: ");
        int dep=sc.nextInt();
        int[] x_bin=decToBin(x1,4);
        int[] y_bin=decToBin(y1,4);
        int[] x_stoch=comparator(rng1,x_bin);
        int[] y_stoch=comparator(rng2,y_bin);
        displaySN(x_stoch);
        displaySN(y_stoch);
        System.out.println("Prev SCC="+scc(x_stoch,y_stoch,16));
        int y_new[]=desynch(x_stoch,y_stoch,dep);
        System.out.print("After desynchronizer: ");
        displaySN(y_new);
        System.out.println("Desynchronizer SCC="+scc(x_stoch,y_new,16));
        int y_new1[]=synch(x_stoch,y_stoch,dep);
        System.out.print("After synchronizer: ");
        displaySN(y_new1);
        System.out.println("Synchronizer SCC="+scc(x_stoch,y_new1,16));
        
        System.out.print("Y/N?");
        char ch=sc.next().charAt(0);
        if (ch=='n'||ch=='N') break;
        }
    }
    public static void main(String args[]){
        Scanner sc=new Scanner(System.in);
        int[] seed1={0,0,0,1};
        int[][] rng1=lfsr(seed1,2,3);
        //for (int x:rng1) System.out.print(x+" ");
        //return;
        int[] seed2={0,0,0,1};
        int[][] rng2=lfsr(seed2,0,3);
        int[] seed3={0,0,0,1};
        int[][] rng3=lfsr(seed3,0,3);
        int px[]={1,0,0,0};
        int py[]={1,0,1,0};
        int sccmag[]={0,0,0,0};
        int sccsign=0;
        /*for (int[] arr:rng2){
            for (int z:arr) System.out.print(z+" ");
            System.out.println();
        }*/
        
        int[] x=comparator(rng1,px); //rng1-->rng2
        System.out.print("Original x:");
        displaySN(x);
        
        int[] y_original=comparator(rng2,py);
        System.out.print("Original y:");
        displaySN(y_original);
        
        int[] sccmagsn=comparator(rng3,sccmag);
        int[] sn01=comparator(rng1,py);
        int[] sn00=comparator(rng2,py);
        int[] sn11=comparator(invert(rng1),py);
        int[] y=new int[16];
        rng1=invert(rng1);
        for (int i=0;i<16;i++){
            if (sccsign==0){
                if (sccmagsn[i]==0) y[i]=sn00[i];
                else y[i]=sn01[i];
            }
            else{
                if (sccmagsn[i]==0) y[i]=sn00[i];
                else y[i]=sn11[i];
            }
        }
        System.out.print("Correlated y: ");
        displaySN(y);
        System.out.println("SCC="+scc(x,y,16));
        
        /*int[] sn1={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};
        int[] sn2={1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0};
        System.out.println("SCC="+scc(sn1,sn2,16));*/
        
        int[] y_zeroSCC=zeroSCC(x,y_original,16);
        System.out.print("Zero SCC y:");
        displaySN(y_zeroSCC);
        System.out.println("SCC="+scc(x,y_zeroSCC,16));
        
        //new code
        
        while(true){
            System.out.println("\n\n***********************************");

            System.out.print("Enter x in decimal: ");
            int x1=sc.nextInt();
            System.out.print("Enter y in decimal: ");
            int y1=sc.nextInt();
            System.out.print("Enter scc: ");
            double s=sc.nextDouble();
            
            int[] x_bin=decToBin(x1,4);
            int[] y_bin=decToBin(y1,4);
            int[] x_stoch=comparator(rng1,x_bin);
            int[] y_stoch=comparator(rng2,y_bin);
            
            
            System.out.println("Required scc="+s);
            System.out.print("Original x:   ");
            displaySN(x_stoch);
            System.out.print("Original y:   ");
            displaySN(y_stoch);
            System.out.println("Original scc: "+scc(x_stoch,y_stoch,16));
            int[] y_anySCC=anySCC(x_stoch,y_stoch,16,s);
            System.out.print("Correlated y:  ");
            displaySN(y_anySCC);
            System.out.println("New y SCC="+scc(x_stoch,y_anySCC,16));
            
            System.out.print("Again?(Y/N): ");
            char query=sc.next().charAt(0);
            if (query=='n'||query=='N') break;
            
        }
    }
}
