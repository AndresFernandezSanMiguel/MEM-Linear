/*============================================================================
 SAOAS - Structural Analysis of Atomic Systems
 Molecular Element Method (MEM): a FEM-based formulation for the linear
 elastic analysis of atomic structures (static response and normal-mode /
 vibrational analysis) from a classical force field.

 If you use this code in your research, please cite:

   [1] A. Fernandez-San Miguel, I. Couceiro, L. Ramirez.
       "The Molecular Element Method (MEM): A FEM-based Formulation for
       Linear and Non-Linear Molecular Elasticity."
       https://ruc.udc.es/entities/publication/d0df64f6-4666-431f-9fbd-0e0f903c5b06

   [2] A. Fernandez-San Miguel, I. Couceiro, L. Ramirez, F. Navarrina.
       "A first order FEM-based formulation for the analysis of molecular
       structures with bonded interactions."
       Engineering with Computers (2024). https://doi.org/10.1007/s00366-024-02085-w

 License: PolyForm Noncommercial 1.0.0 - see LICENSE file.
 Repository: https://github.com/<your-username>/MEM-saoas
==============================================================================*/

#include<iostream>
#include<iomanip>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include<cmath>
#include<Eigen/Core>
#include<Eigen/Geometry>
#include<Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include<fstream>
#include<sstream>
#include<vector>
#include<Spectra/SymEigsSolver.h>
#include<Spectra/SymGEigsSolver.h>
#include<Spectra/MatOp/SparseSymMatProd.h>
#include<Spectra/MatOp/SparseCholesky.h>
#include<Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace Spectra;

/*Estructuras de Conectividad__________________*/

struct Barra{ int i,j; double kl; };
struct Angulo{ int i,j,k; double ko; };
struct Diedro{ int i,j,k,l; double kp; };

/*Clase Operadora original (K definida positiva, con restricciones)_________________________*/
// Recibe un solver ya factorizado para reutilizar la factorizacion LDLT
// (evita factorizar K dos veces cuando tambien se resuelve el sistema estatico)
class ShiftInvertOp {
public:
    using Scalar = double;
private:
    const SparseMatrix<double>& m_M;
    const SimplicialLDLT<SparseMatrix<double>>& m_solver;
public:
    ShiftInvertOp(const SimplicialLDLT<SparseMatrix<double>>& solver,
                  const SparseMatrix<double>& M)
        : m_M(M), m_solver(solver) {}
    int rows() const { return m_M.rows(); }
    int cols() const { return m_M.cols(); }
    void perform_op(const double* x_in, double* y_out) const {
        int n = rows();
        Map<const VectorXd> x(x_in, n);
        VectorXd y = m_solver.solve(m_M * x);
        Map<VectorXd>(y_out, n) = y;
    }
};

/*Menu__________________________________________*/

void lectura(int &nnodo,int &neles,int &nelan,int &neldi,string &modo);

void dim(int nnodo,int neles,int nelan,int neldi,int &nres,
          const string &modo,
          Matrix<double,4,Dynamic> &mcoord,
          vector<Barra> &v_bar,vector<Angulo> &v_ang,vector<Diedro> &v_die,
          MatrixXi &mres);

void ensames(const vector<Barra> &v_bar,MatrixXi &mres,
             vector<Triplet<double>> &kg_triplets,const Matrix<double,4,Dynamic> &mcoord);

void ensaman(const vector<Angulo> &v_ang,MatrixXi &mres,
             vector<Triplet<double>> &kg_triplets,const Matrix<double,4,Dynamic> &mcoord);

void ensamdi(const vector<Diedro> &v_die,MatrixXi &mres,
             vector<Triplet<double>> &kg_triplets,const Matrix<double,4,Dynamic> &mcoord);

void ensamgen(const vector<int> &nodos,MatrixXi &mres,
               MatrixXd &kel,vector<Triplet<double>> &kg_triplets);

void calcul(int nnodo,MatrixXi &mres,
             const Matrix<double,4,Dynamic> &mcoord,SparseMatrix<double> &kg,
             SparseMatrix<double> &masm,const string &modo);

void estatico(int nnodo,MatrixXi &mres,
              const SimplicialLDLT<SparseMatrix<double>>& solver,
              const SparseMatrix<double>& kg,
              const Matrix<double,4,Dynamic>& mcoord);

void normal(int nnodo,double lamlim,const VectorXd &evalues,
            MatrixXi &mres,const MatrixXd &evecs,const Matrix<double,4,Dynamic> &mcoord);

/*Programa Principal-------------------------*/

int main(){
    int nnodo{},neles{},nelan{},neldi{},nres{};
    string modo{};

    lectura(nnodo,neles,nelan,neldi,modo);
    
    Matrix<double,4,Dynamic> mcoord(4,nnodo);
    MatrixXi mres(2,3*nnodo);
    
    vector<Barra> v_bar;
    vector<Angulo> v_ang;
    vector<Diedro> v_die;

    dim(nnodo,neles,nelan,neldi,nres,modo,mcoord,v_bar,v_ang,v_die,mres);

    int dim_reducida=3*nnodo-nres;
    SparseMatrix<double> kg(dim_reducida,dim_reducida);
    SparseMatrix<double> masm(dim_reducida,dim_reducida);

    vector<Triplet<double>> kg_triplets;
    kg_triplets.reserve(neles*21+nelan*45+neldi*78);

    ensames(v_bar,mres,kg_triplets,mcoord);
    ensaman(v_ang,mres,kg_triplets,mcoord);
    ensamdi(v_die,mres,kg_triplets,mcoord);

    kg.setFromTriplets(kg_triplets.begin(),kg_triplets.end());
    kg=SparseMatrix<double>(kg.selfadjointView<Upper>());

    calcul(nnodo,mres,mcoord,kg,masm,modo);

    return 0;
}

/*Implementacion de Funciones__________________________________________*/

void lectura(int &nnodo,int &neles,int &nelan,int &neldi,string &modo){
    ifstream le("datos.txt");
    if(!le){ cout<<"Error: datos.txt not found."<<endl; return; }
    string linea;
    auto skip=[&](int n){ for(int i=0;i<n;i++) getline(le,linea); };
    skip(9);
    getline(le,linea);
    if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
    { auto p=linea.find('='); if(p!=string::npos) modo=linea.substr(p+1); }
    skip(11);
    for(int i=0;i<4;i++){
        getline(le,linea);
        if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
        auto p=linea.find('='); if(p==string::npos) continue;
        string clave=linea.substr(0,p); int valor=stoi(linea.substr(p+1));
        if(clave=="nnodo")      nnodo=valor;
        else if(clave=="neles") neles=valor;
        else if(clave=="nelan") nelan=valor;
        else if(clave=="neldi") neldi=valor;
    }
    le.close();
}

void dim(int nnodo,int neles,int nelan,int neldi,int &nres,
          const string &modo,
          Matrix<double,4,Dynamic> &mcoord,
          vector<Barra> &v_bar,vector<Angulo> &v_ang,vector<Diedro> &v_die,
          MatrixXi &mres){

    nres=0;
    for(int i=0;i<3*nnodo;i++){ mres(0,i)=1; mres(1,i)=i; }

    bool con_res=(modo=="RV"||modo=="S");
    double rig{};
    int i_idx{},j_idx{},k_idx{},l_idx{};

    ifstream le("datos.txt");
    if(!le){ cout<<"Error: datos.txt not found."<<endl; return; }
    string linea;
    auto skip=[&](int n){ for(int i=0;i<n;i++) getline(le,linea); };

    // 13 lineas entre neldi= y primer dato atomico
    skip(38);
    for(int i=0;i<nnodo;i++){
        getline(le,linea);
        if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
        istringstream ss(linea);
        int ni; double m,x,y,z;
        ss>>ni>>m>>x>>y>>z;
        mcoord(0,i)=x; mcoord(1,i)=y;
        mcoord(2,i)=z; mcoord(3,i)=m;
        if(con_res){
            int u,v,w; double fx,fy,fz;
            ss>>u>>v>>w>>fx>>fy>>fz;
            int dofs[3]={u,v,w};
            for(int j=0;j<3;j++){
                mres(0,3*i+j)=dofs[j];
                if(dofs[j]==0) mres(1,3*i+j)=-10;
            }
        }
    }
    skip(11);
    for(int i=0;i<neles;i++){
        getline(le,linea);
        if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
        istringstream ss(linea); ss>>i_idx>>j_idx>>rig;
        v_bar.push_back({i_idx,j_idx,rig});
    }
    skip(11);
    for(int i=0;i<nelan;i++){
        getline(le,linea);
        if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
        istringstream ss(linea); ss>>i_idx>>j_idx>>k_idx>>rig;
        v_ang.push_back({i_idx,j_idx,k_idx,rig});
    }
    skip(11);
    for(int i=0;i<neldi;i++){
        getline(le,linea);
        if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
        istringstream ss(linea); ss>>i_idx>>j_idx>>k_idx>>l_idx>>rig;
        v_die.push_back({i_idx,j_idx,k_idx,l_idx,rig});
    }
    le.close();

    if(con_res){
        nres=0;
        for(int i=0;i<3*nnodo;i++){
            if(mres(0,i)==0){ mres(1,i)=-10; ++nres; }
            else mres(1,i)=i-nres;
        }
    }
}

void ensames(const vector<Barra> &v_bar,MatrixXi &mres,
             vector<Triplet<double>> &kg_triplets,const Matrix<double,4,Dynamic> &mcoord){
    MatrixXd kel(6,6);
    for(const auto &b : v_bar){
        Vector3d rij=mcoord.col(b.j).head<3>()-mcoord.col(b.i).head<3>();
        rij.normalize();
        VectorXd grad(6);
        grad<<-rij,rij;
        kel=(b.kl*grad)*grad.transpose();
        ensamgen({b.i,b.j},mres,kel,kg_triplets);
    }
}

void ensaman(const vector<Angulo> &v_ang,MatrixXi &mres,
             vector<Triplet<double>> &kg_triplets,const Matrix<double,4,Dynamic> &mcoord){
    MatrixXd kel(9,9);
    for(const auto &a : v_ang){
        Vector3d rji=mcoord.col(a.i).head<3>()-mcoord.col(a.j).head<3>();
        Vector3d rjk=mcoord.col(a.k).head<3>()-mcoord.col(a.j).head<3>();
        Vector3d vn=rji.cross(rjk);
        vn.normalize();
        double lrji=rji.norm(); double lrjk=rjk.norm();
        Vector3d g1=(rji/lrji).cross(vn)/lrji;
        Vector3d g3=vn.cross(rjk/lrjk)/lrjk;
        VectorXd grad(9);
        grad<<g1,-(g1+g3),g3;
        kel=(a.ko*grad)*grad.transpose();
        ensamgen({a.i,a.j,a.k},mres,kel,kg_triplets);
    }
}

void ensamdi(const vector<Diedro> &v_die,MatrixXi &mres,
             vector<Triplet<double>> &kg_triplets,const Matrix<double,4,Dynamic> &mcoord){
    MatrixXd kel(12,12);
    for(const auto &d : v_die){
        Vector3d r21=mcoord.col(d.i).head<3>()-mcoord.col(d.j).head<3>();
        Vector3d r32=mcoord.col(d.j).head<3>()-mcoord.col(d.k).head<3>();
        Vector3d r34=mcoord.col(d.l).head<3>()-mcoord.col(d.k).head<3>();
        Vector3d vn1=r21.cross(r32); Vector3d vn2=r34.cross(r32);
        double lr21=r21.norm(); double lr32=r32.norm(); double lr34=r34.norm();
        double lvn1=vn1.norm(); double lvn2=vn2.norm();
        double c1=-(1.0+lr21/lr32*((r21/lr21).dot(r32/lr32)));
        double c2=-lr34/lr32*((r34/lr34).dot(r32/lr32));
        Vector3d g1=-lr32/lvn1*(vn1/lvn1);
        Vector3d g4=lr32/lvn2*(vn2/lvn2);
        Vector3d g2=c1*g1+c2*g4; Vector3d g3=-(g1+g2+g4);
        VectorXd grad(12);
        grad<<g1,g2,g3,g4;
        kel=(d.kp*grad)*grad.transpose();
        ensamgen({d.i,d.j,d.k,d.l},mres,kel,kg_triplets);
    }
}

void ensamgen(const vector<int> &nodos,MatrixXi &mres,
               MatrixXd &kel,vector<Triplet<double>> &kg_triplets){
    int nt=nodos.size();
    for(int i=0;i<nt;i++){
        for(int j=0;j<nt;j++){
            for(int k=0;k<3;k++){
                for(int l=0;l<3;l++){
                    int f_g=3*nodos[i]+k;
                    int c_g=3*nodos[j]+l;
                    if(mres(0,f_g)!=0 && mres(0,c_g)!=0){
                        int f_map=mres(1,f_g);
                        int c_map=mres(1,c_g);
                        if(f_map<=c_map) kg_triplets.push_back(Triplet<double>(f_map,c_map,kel(3*i+k,3*j+l)));
                    }
                }
            }
        }
    }
}

void calcul(int nnodo,MatrixXi &mres,
             const Matrix<double,4,Dynamic> &mcoord,SparseMatrix<double> &kg,
             SparseMatrix<double> &masm,const string &modo){
    int aumax{}; const double cluz{0.000299792458}; double lamlim{};

    // Ensamblar matriz de masas diagonal
    masm.reserve(VectorXi::Constant(masm.rows(),1));
    for(int i=0;i<nnodo;i++){
        for(int j=0;j<3;j++){
            if(mres(0,3*i+j)!=0){
                int idx=mres(1,3*i+j);
                masm.insert(idx,idx)=mcoord(3,i);
            }
        }
    }
    masm.makeCompressed();

    bool molecula_libre = (kg.rows() == 3*nnodo);

    // Mode S: static analysis only, no eigenvalues
    if(modo=="S"){
        SimplicialLDLT<SparseMatrix<double>> solver_k;
        solver_k.compute(kg);
        if(solver_k.info() != Success){
            cout<<"Error: LDLT factorization of K failed."<<endl;
            return;
        }
        estatico(nnodo, mres, solver_k, kg, mcoord);
        return;
    }

    cout<<"Enter number of eigenvalues: "; cin>>aumax;
    cout<<"Enter lambda threshold: "; cin>>lamlim;

    VectorXd evalues;
    MatrixXd evecs;

    if(molecula_libre){
        // --- Molecula libre: problema generalizado K v = lambda M v (modo Cholesky) ---
        // K es semidefinida positiva (6 modos rigidos con lambda=0).
        // SparseSymMatProd: y = K*x, sin factorizacion (K puede ser semidefinida).
        // SparseCholesky:   factoriza M, siempre definida positiva con masas > 0.
        // Spectra aplica M^{-1/2} K M^{-1/2} on-the-fly sin materializar A.
        // LargestMagn: los k modos de mayor autovalor quedan seleccionados;
        // los modos rigidos (lambda=0) quedan excluidos automaticamente.
        cout<<"Free molecule: generalized eigenproblem K v=lambda M v (Cholesky mode)."<<endl;

        int n   = kg.rows();
        int nev = min(aumax, n-1);
        int ncv = min(max(2*nev, nev+10), n-1);
        if(nev < 1){ cout<<"ERROR: System too small."<<endl; return; }

        SparseSymMatProd<double> opK(kg);
        SparseCholesky<double>   opM(masm);

        SymGEigsSolver<SparseSymMatProd<double>,
                       SparseCholesky<double>,
                       GEigsMode::Cholesky> geigs(opK, opM, nev, ncv);
        geigs.init();
        int nconv = geigs.compute(SortRule::LargestMagn, 500, 1e-14);

        if(geigs.info() == CompInfo::Successful){
            evalues = geigs.eigenvalues();
            evecs   = geigs.eigenvectors();
            cout<<"Converged eigenvalues: "<<nconv<<endl;
        } else {
            cout<<"ERROR: SymGEigsSolver did not converge"
                <<" (info="<<(int)geigs.info()<<", converged: "
                <<nconv<<" of "<<nev<<")."<<endl;
            cout<<"Try requesting fewer eigenvalues."<<endl;
            return;
        }

    } else {
        // --- Con restricciones: K definida positiva, Spectra sparse ---
        

        // Factorizamos K una sola vez: sirve tanto para el sistema estatico
        // como para el operador shift-invert de vibraciones
        SimplicialLDLT<SparseMatrix<double>> solver_k;
        solver_k.compute(kg);
        if(solver_k.info() != Success){
            cout<<"ERROR: LDLT factorization of K failed."<<endl;
            return;
        }

        ShiftInvertOp op(solver_k, masm);
        int n = op.rows();
        int nev = min(aumax, n-1);
        int ncv = min(max(2*nev, nev+10), n-1);
        if(nev < 1){ cout<<"ERROR: System too small."<<endl; return; }

        SymEigsSolver<ShiftInvertOp> geigs(op, nev, ncv);
        geigs.init();
        int nconv = geigs.compute(SortRule::LargestMagn, 500, 1e-14);

        if(geigs.info() == CompInfo::Successful){
            VectorXd ritz = geigs.eigenvalues();
            evalues.resize(ritz.size());
            for(int i=0;i<ritz.size();i++) evalues(i) = 1.0/ritz(i);
            evecs = geigs.eigenvectors();
            cout<<"Converged eigenvalues: "<<nconv<<endl;
        } else {
            cout<<"ERROR: Spectra solver failed (info="<<(int)geigs.info()<<")"<<endl;
            return;
        }
    }

    // Write autoval.txt
    ofstream feig("autoval.txt");
    if(feig){
        feig<<"________________________________________________________"<<endl;
        feig<<"vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"<<endl;
        feig<<"i<--Number of the eigenvalue"<<endl;
        feig<<"lambda_i<--Value of the i-th eigenvalue"<<endl;
        feig<<"--------------------------------------------------------"<<endl;
        feig<<"i lambda_i"<<endl;
        int cont_eig=0;
        for(int i=0;i<evalues.size();i++){
            if(evalues(i)>lamlim){
                ++cont_eig;
                feig<<cont_eig<<" "<<scientific<<setprecision(12)<<evalues(i)<<endl;
            }
        }
    }
    feig.close();
    normal(nnodo,lamlim,evalues,mres,evecs,mcoord);
}

void estatico(int nnodo,MatrixXi &mres,
              const SimplicialLDLT<SparseMatrix<double>>& solver,
              const SparseMatrix<double>& kg,
              const Matrix<double,4,Dynamic>& mcoord){
    // Read forces from datos.txt (columns f_xi f_yi f_zi)
    VectorXd fuerza_global = VectorXd::Zero(3*nnodo);
    ifstream lef("datos.txt");
    if(!lef){ cout<<"Error: cannot reread datos.txt"<<endl; return; }
    string linea;
    auto skip=[&](int n){ for(int i=0;i<n;i++) getline(lef,linea); };
    skip(38);
    for(int i=0;i<nnodo;i++){
        getline(lef,linea);
        if(!linea.empty()&&linea.back()=='\r') linea.pop_back();
        istringstream ss(linea);
        int ni,u,v,w; double m,x,y,z,fx,fy,fz;
        ss>>ni>>m>>x>>y>>z>>u>>v>>w>>fx>>fy>>fz;
        fuerza_global(3*i)  =fx;
        fuerza_global(3*i+1)=fy;
        fuerza_global(3*i+2)=fz;
    }
    lef.close();

    // Construir vector de fuerzas reducido (solo DOFs activos)
    int n = kg.rows();
    VectorXd f(n);
    for(int i=0;i<3*nnodo;i++){
        if(mres(0,i)!=0){
            int idx=mres(1,i);
            f(idx)=fuerza_global(i);
        }
    }

    // Resolver K*u=f con la factorizacion LDLT ya disponible
    VectorXd u_red = solver.solve(f);
    if(solver.info() != Success){
        cout<<"Error al resolver el sistema estatico."<<endl;
        return;
    }

    // Expandir resultado: insertar ceros en DOFs restringidos
    VectorXd u_global = VectorXd::Zero(3*nnodo);
    for(int i=0;i<3*nnodo;i++){
        if(mres(0,i)!=0) u_global(i)=u_red(mres(1,i));
    }

    // Write displacements to disp.txt
    ofstream escr("disp.txt");
    if(escr){
        escr<<"________________________________________________________"<<endl;
        escr<<"vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"<<endl;
        escr<<"N<--Total number of nuclei"<<endl;
        escr<<"u<--Global displacement vector (dim=3*N)"<<endl;
        escr<<"i<--Index of u"<<endl;
        escr<<"N_i<--Global numeration of the nuclei"<<endl;
        escr<<"u_i<--i-th component of u"<<endl;
        escr<<"--------------------------------------------------------"<<endl;
        escr<<"i N_i u_i"<<endl;
        for(int i=0;i<nnodo;i++){
            for(int j=0;j<3;j++){
                escr<<j<<" "<<i<<" "<<scientific<<setprecision(12)<<u_global(3*i+j)<<endl;
            }
        }
        cout<<"Displacements written to disp.txt"<<endl;
    }
    escr.close();

    // Write deformation.pdb: multi-model PDB trajectory from original to deformed
    // B-factor = displacement magnitude |u_i| for color mapping in VMD

    // Precompute displacement magnitudes, normalized to [0,100] for B-factor
    // VMD colors by B-factor; [0,100] gives full color range
    VectorXd bmag(nnodo);
    for(int i=0;i<nnodo;i++){
        double ux=u_global(3*i), uy=u_global(3*i+1), uz=u_global(3*i+2);
        bmag(i)=sqrt(ux*ux+uy*uy+uz*uz);
    }
    double bmax=bmag.maxCoeff();
    if(bmax>0) bmag = (bmag/bmax)*99.99;

    // Loop: allow user to regenerate PDB with different scale/frames without recomputing
    char again='y';
    while(again=='y'||again=='Y'){
        double escal_def{};
        int nframes{};
        cout<<"Enter displacement scale factor for deformation file: "; cin>>escal_def;
        cout<<"Enter number of frames (original->deformed): "; cin>>nframes;
        if(nframes<2) nframes=2;

        ofstream fdef("deformation.pdb");
        if(fdef){
            fdef<<fixed<<setprecision(3);
            for(int t=0;t<nframes;t++){
                double alpha=(nframes>1)?(double)t/(nframes-1):0.0;
                fdef<<"MODEL     "<<setw(4)<<t+1<<endl;
                for(int i=0;i<nnodo;i++){
                    double x=mcoord(0,i)+escal_def*alpha*u_global(3*i);
                    double y=mcoord(1,i)+escal_def*alpha*u_global(3*i+1);
                    double z=mcoord(2,i)+escal_def*alpha*u_global(3*i+2);
                    // Strict PDB column format
                    char pdbline[81];
                    snprintf(pdbline,81,
                        "ATOM  %5d  CA  MOL A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          C ",
                        i+1, 1, x, y, z, 1.00, bmag(i));
                    fdef<<pdbline<<endl;
                }
                fdef<<"ENDMDL"<<endl;
            }
            cout<<"Deformation trajectory written to deformation.pdb ("<<nframes<<" frames)"<<endl;
        }
        fdef.close();
        cout<<"Generate again with different scale? (y/n): "; cin>>again;
    }
}

void normal(int nnodo,double lamlim,const VectorXd &evalues,
            MatrixXi &mres,const MatrixXd &evecs,const Matrix<double,4,Dynamic> &mcoord){
    int cont{}; double escal{};
    cout<<"Enter scale factor: "; cin>>escal;
    ofstream escr("normal.nmd");
    if(escr){
        escr<<"title xxxx \nnames ";
        for(int i=0;i<nnodo;i++) escr<<"C"<<i+1<<" ";
        escr<<"\nresnames ";
        for(int i=0;i<nnodo;i++) escr<<"X ";
        escr<<"\nresnums ";
        for(int i=0;i<nnodo;i++) escr<<"1 ";
        escr<<"\ncoordinates ";
        for(int i=0;i<nnodo;i++){
            for(int j=0;j<3;j++) escr<<setprecision(12)<<escal*mcoord(j,i)<<" ";
        }
        escr<<"\n";
        for(int i=0;i<evalues.size();i++){ 
            if(evalues(i)>lamlim){
                int contres=0;
                escr<<"mode "<<++cont<<" ";
                for(int j=0;j<3*nnodo;j++){
                    if(mres(0,j)==0){ ++contres; escr<<0<<" "; }
                    else escr<<setprecision(12)<<evecs(j-contres,i)<<" ";
                }
                escr<<"\n";
            }
        }
    }
    escr.close();
}