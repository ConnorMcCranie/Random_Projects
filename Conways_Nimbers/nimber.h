// nimber - GF[2^2^k] as described by J.H. Conway
// David Eppstein, Columbia University, 19 Mar 1988

typedef unsigned long proto_nimber;

class nimber {
    proto_nimber val;

public:
    nimber(proto_nimber x) { val = x; }
    proto_nimber value() { return val; }

    int operator== (nimber that) { return this->val == that.val; }
    int operator!= (nimber that) { return this->val != that.val; }

    nimber square();
    nimber sqrt();
    nimber inverse();

    nimber operator* (nimber);
    nimber operator+ (nimber that) { return nimber(this->val ^ that.val); }
    nimber operator/ (nimber that) { return (*this) * that.inverse(); }
};
