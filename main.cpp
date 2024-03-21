#include <pybind11/pybind11.h>
#include <bits/stdc++.h>
using namespace std;
using ld = long double;
namespace py = pybind11;
class Polynomial {
private:
    vector<ld> arr;
    void addT(ld val, int pow){
        while(arr.size()<=pow) arr.push_back(0);
        arr[pow]+=val;
    }
    vector<complex<ld>> fft(vector<complex<ld>> curr, bool inv){
        int n = curr.size();
        if(n == 1) return {complex<ld>(curr[0])};
        vector<complex<ld>> even, odd;
        for(int i = 0; i<n; i++)
            (i&1?odd:even).push_back(curr[i]);
        vector<complex<ld>> efft = fft(even,inv);
        vector<complex<ld>> offt = fft(odd,inv);
        complex<ld> w = polar((ld)1.0,(ld)(inv?-2:2)*3.141592653589793/n);
        vector<complex<ld>> ans(n);
        complex<ld> pow = complex<ld>(1);
        for(int i = 0; i<n>>1; i++){
            ans[i] = efft[i]+pow*offt[i];
            ans[i+(n>>1)] = efft[i]-pow*offt[i];
            pow*=w;
        }
        return ans;
    }
public:
    Polynomial(){}
    Polynomial(const Polynomial& o) {
        for(const ld& a: o.arr) arr.push_back(a);
    }
    Polynomial(py::list list){
        for(auto& a: list) arr.push_back(a.cast<double>());
    }
    py::float_ operator[](py::object val){
        ld d = val.cast<ld>();
        ld x = 1, ans = 0;
        for(ld& a: arr){
            ans+=x*a;
            x*=d;
        }
        py::float_ p((double)ans);
        return p;
    }
    py::float_ getTerm(int pow){
        if(pow>=arr.size()) return py::float_(0.0);
        return py::float_((double)arr[pow]);
    }
    void setTerm(py::float_ val, int pow){
        assert(pow>=0);
        while(arr.size()<=pow) arr.push_back(0);
        arr[pow] = (ld)val.cast<double>();
    }
    void addTerm(py::float_ val, int pow){
        assert(pow>=0);
        while(arr.size()<=pow) arr.push_back(0);
        arr[pow]+=val.cast<double>();
    }
    Polynomial operator+(Polynomial& other) {
        Polynomial ans;
        for(int i = 0; i<arr.size(); i++) ans.addT(arr[i],i);
        for(int i = 0; i<other.arr.size(); i++) other.addT(other.arr[i],i);
        return ans;
    }
    Polynomial operator*(Polynomial& other) {
        int sz = arr.size()+other.arr.size();
        if(sz >= 50){
            vector<complex<ld>> x, y;
            for(ld a: arr) x.push_back(a);
            for(ld a: other.arr) y.push_back(a);
            while(__builtin_popcount(sz) != 1) sz++;
            while(x.size()<sz) x.push_back(complex<ld>(0));
            while(y.size()<sz) y.push_back(complex<ld>(0));
            vector<complex<ld>> a = fft(x,0);
            vector<complex<ld>> b = fft(y,0);
            vector<complex<ld>> mul(a.size());
            for(int i = 0; i<a.size(); i++) mul[i] = a[i]*b[i];
            vector<complex<ld>> ans = fft(mul,1);
            Polynomial ret;
            for(int i = 0; i<ans.size(); i++) ret.arr.push_back(ans[i].real()/ans.size());
            return ret;
        }else{
            Polynomial ans;
            for(int i = 0; i<arr.size(); i++)
                for(int j = 0; j<other.arr.size(); j++)
                    ans.addT(arr[i]*other.arr[i],i+j);
            return ans;
        }
    }
    Polynomial operator*(ld other){
        Polynomial ans;
        for(int i = 0; i<arr.size(); i++)
            ans.addT(arr[i]*other,i);
        return ans;
    }
    Polynomial operator-(Polynomial& other) {
        Polynomial ans;
        for(int i = 0; i<arr.size(); i++) ans.addT(arr[i],i);
        for(int i = 0; i<other.arr.size(); i++) other.addT(-other.arr[i],i);
        return ans;
    }
    Polynomial operator/(ld other){
        assert(other != 0);
        Polynomial ans;
        for(int i = 0; i<arr.size(); i++)
            ans.addT(arr[i]/other,i);
        return ans;
    }
    Polynomial& operator+=(Polynomial& other) {
        for(int i = 0; i<other.arr.size(); i++) addT(other.arr[i],i);
        return *this;
    }
    Polynomial& operator-=(Polynomial& other) {
        for(int i = 0; i<other.arr.size(); i++) addT(-other.arr[i],i);
        return *this;
    }
    Polynomial& operator*=(ld other){
        for(int i = 0; i<arr.size(); i++)
            arr[i]*=other;
        return *this;
    }
    int degree() const {
        return arr.size()-1;
    }
    pair<Polynomial,Polynomial> div(Polynomial& other) const {
        Polynomial curr, mul;
        for(ld i: arr) curr.arr.push_back(i);
        int m = other.degree();
        while(m >= 0 && abs(other.arr[m])<1e-8) m--;
        assert(m != -1);
        for(int i = degree(); i>=m; i--){
            ld div = curr.arr[i]/other.arr[m];
            for(int j = m; j>=0; j--) curr.arr[i+j-m]-=div*other.arr[j];
            mul.addT(div,i-m);
        }
        return make_pair(mul,curr);
    }
};
PYBIND11_MODULE(polynomial, handle){
    handle.doc() = "Test docs";
    py::class_<Polynomial>(handle,"Polynomial")
        .def(py::init<py::list>())
        .def(py::init<>())
        .def(py::init<Polynomial>())
        .def("__getitem__",&Polynomial::operator[])
        .def("__add__",&Polynomial::operator+)
        .def("__mul__",[](Polynomial& a, Polynomial& b){return a*b;})
        .def("__mul__",[](Polynomial& a, double b){return a*b;})
        .def("__sub__",&Polynomial::operator-)
        .def("__mod__", [](Polynomial& a, Polynomial& b){return move(a.div(b).second);})
        .def("__floordiv__", [](Polynomial& a, Polynomial& b){return move(a.div(b).first);})
        .def("__iadd__",&Polynomial::operator+=)
        .def("__isub__",&Polynomial::operator-=)
        .def("__imul__",&Polynomial::operator*=)
        .def("__imul__",[](Polynomial& a, Polynomial b){return a = a*b;})
        .def("__ifloordiv__",[](Polynomial& a, Polynomial b){return a = a.div(b).first;})
        .def("__imod__",[](Polynomial& a, Polynomial b){return a = a.div(b).second;})
        .def("__div__",&Polynomial::div)
        .def("degree",&Polynomial::degree)
        .def("getTerm",&Polynomial::getTerm)
        .def("setTerm",&Polynomial::setTerm)
        .def("addTerm",&Polynomial::addTerm)
        ;
}