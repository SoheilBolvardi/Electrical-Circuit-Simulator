#include <bits/stdc++.h>
#include <fstream>
#include <sys/stat.h>
#include <direct.h>
#include <unistd.h>
#include <windows.h>
#include <dirent.h>


using namespace std;

vector<double> gaussianElimination(const vector<vector<double>> &A_in, const vector<double> &b_in) {
    int n = b_in.size();
    vector<vector<double>> A = A_in;
    vector<double> b = b_in;
    vector<double> x(n, 0.0);

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int row = i + 1; row < n; ++row) {
            if (abs(A[row][i]) > abs(A[pivot][i])) {
                pivot = row;
            }
        }
        if (abs(A[pivot][i]) < 1e-12) {
            throw runtime_error("Matrix is singular or nearly singular");
        }

        if (pivot != i) {
            swap(A[i], A[pivot]);
            swap(b[i], b[pivot]);
        }

        for (int row = i + 1; row < n; ++row) {
            double factor = A[row][i] / A[i][i];
            for (int col = i; col < n; ++col) {
                A[row][col] -= factor * A[i][col];
            }
            b[row] -= factor * b[i];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int col = i + 1; col < n; ++col) {
            sum -= A[i][col] * x[col];
        }
        x[i] = sum / A[i][i];
    }

    return x;
}

double clean(double x) {
    if (abs(x) < 1e-5) return 0.0;
    return x;
}

class Node;

class Element;

class Node {
private:
    bool isGround;
    string name;
    double voltage;
    vector<Element *> connected_elements;
public:
    Node(string name_) : name(name_), isGround(false) {}

    vector<Element *> getConnectedElements() { return connected_elements; }

    string getName() { return name; }

    void addConnectedElement(Element *e) {
        connected_elements.push_back(e);
    }

    void removeConnectedElement(Element *e) {
        auto x = find(connected_elements.begin(), connected_elements.end(), e);
        if (x != connected_elements.end())
            connected_elements.erase(x);
    }

    void setName(string new_name) {
        name = new_name;
    }

    void setVoltage(double v) { voltage = v; }

    double getVoltage() { return voltage; }

    void setGround() { isGround = true; }

    bool isGroundd() { return isGround; }
};

class Element {
protected:
    string name;
    double value;
    Node *node1;
    Node *node2;
    double current = 0;
    double currenttime = 0;
public:
    Element(string name_, double value_, Node *n1, Node *n2) :
            name(name_), value(value_), node1(n1), node2(n2) {}

    string getName() {
        return name;
    }

    void setValue(double value_) {
        value = value_;
    }

    virtual string getType() = 0;

    virtual string getCategory() { return "Generic"; }

    virtual void stamp(vector<vector<double>> &A,
                       vector<double> &b,
                       map<string, int> &nodeToIndex,
                       int &voltageIndex,
                       map<Element *, int> &currentIndexmap) = 0;

    virtual double getValue() { return value; }

    Node *getFirstNode() { return node1; }

    Node *getSecondNode() { return node2; }

    virtual double getCurrent() { return current; }

    virtual double getVoltage() {}

    virtual void setCurrent(double i) { current = i; }

    virtual void updateTime(double t) { currenttime = t; }

    virtual ~Element() = default;
};

class Capacitor : public Element {
private:
    double volt_last=0;
public:
    Capacitor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Capacitor"; }

    double getCurrent() override {
        return 0.0;
    }

    double getVoltage() override {

    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        double dt = 1e-3;
        double g = getValue() / dt;
        double v_prev = volt_last;
        double i_eq = g * v_prev;

        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();
        int i = hasN1 ? nodeToIndex[node1->getName()] : -1;
        int j = hasN2 ? nodeToIndex[node2->getName()] : -1;

        if (hasN1) {
            A[i][i] += g;
            if (hasN2) {
                A[i][j] -= g;
                A[j][i] -= g;
                A[j][j] += g;
            }
            b[i] += i_eq;
            if (hasN2) b[j] -= i_eq;
        } else if (hasN2) {
            A[j][j] += g;
            b[j] -= i_eq;
        }
    }



    void updateVoltage(double voltage){
        volt_last=voltage;
    }

    void DCstamp(vector<vector<double>> &A,
                 vector<double> &b,
                 map<string, int> &nodeToIndex,
                 int &voltageIndex,
                 map<Element *, int> &currentIndexmap){

    }
    double getVoltLast() {return volt_last;}
};

class Inductor : public Element {
private:
    double curr_last=0;
    double volt_last=0;
public:
    Inductor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Inductor"; }

    double getCurrent() override {
        double dv = node1->getVoltage() - node2->getVoltage();
        return 1e9 * dv;
    }

    double getVoltage() override {

    }

    string getCategory() override { return "Inductor"; }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        double dt = 1e-3;
        double g = value / dt;
        int vIdx = currentIndexmap[this];

        string n1 = node1->getName();
        string n2 = node2->getName();
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();

        int i = hasN1 ? nodeToIndex[n1] : -1;
        int j = hasN2 ? nodeToIndex[n2] : -1;

        if (hasN1 && hasN2) {
            A[i][vIdx] += 1;
            A[j][vIdx] -= 1;
            A[vIdx][i] += 1;
            A[vIdx][j] -= 1;
        } else if (hasN1) {
            A[i][vIdx] += 1;
            A[vIdx][i] += 1;
        } else if (hasN2) {
            A[j][vIdx] -= 1;
            A[vIdx][j] -= 1;
        }

        A[vIdx][vIdx] -= g;
        b[vIdx] -= curr_last * g;
    }

    void DCstamp(vector<vector<double>> &A,
                 vector<double> &b,
                 map<string, int> &nodeToIndex,
                 int &voltageIndex,
                 map<Element *, int> &currentIndexmap){
        int i = nodeToIndex[getFirstNode()->getName()];
        int j = nodeToIndex[getSecondNode()->getName()];
        string n1= node1->getName();
        string n2= node2->getName();
        double g = 1e9;
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();
        if (hasN1 && hasN2) {
            int ip = nodeToIndex[n1];
            int jp = nodeToIndex[n2];
            A[ip][ip] += g;
            A[jp][jp] += g;
            A[ip][jp] -= g;
            A[jp][ip] -= g;
        } else if (hasN1) {
            int ip = nodeToIndex[n1];
            A[ip][ip] += g;
        } else if (hasN2) {
            int jp = nodeToIndex[n2];
            A[jp][jp] += g;
        }
    }
    void updateCurrent(double i){
        curr_last= i;
    }

    void updateCurrentFromVoltage(double vL, double dt) {
        curr_last += (dt / value) * vL;
    }

    double getCurrLast() {return curr_last;}
};


class Circuit {
protected:
    double deltat=1e-6;
    vector<Element *> elements;
    vector<Node *> nodes;
    map<string, Node *> node_access;
    Node *gnd = nullptr;
    vector<vector<double>> A;
    vector<double> b;
    string Analysistype;
public:
    map<Element *, int> currentIndexmap;

    vector<Element *> &getElements() { return elements; }

    map<string, Node *> getNodeAccess() { return node_access; }

    map<string, int> nodetoindex;
    vector<Node *> ordernodes;

    void setDeltat(double dt){deltat=dt;}
    double getDeltat() {return deltat;}

    void setType(string type){Analysistype=type;}


    Node *getCreateNode(string name) {
        if (!node_access.count(name)) {
            Node *newnode = new Node(name);
            nodes.push_back(newnode);
            node_access[name] = newnode;
        }
        return node_access[name];
    }

    void addElement(Element *element) {
        elements.push_back(element);
    }

    Node *getGround() { return gnd; }

    vector<Node *> getNodes() { return nodes; }

    bool setGround(string name) {
        if (!node_access.count(name))
            return false;
        gnd = node_access[name];
        return true;
    }

    bool removeGround(string name) {
        if (gnd != node_access[name])
            return false;
        gnd = nullptr;
        return true;
    }

    bool renameNode(string old_name, string new_name, string &err) {
        if (node_access.find(old_name) == node_access.end()) {
            err = "ERROR: Node " + old_name + " does not exist in the circuit\n";
            return false;
        }
        if (node_access.find(new_name) != node_access.end()) {
            err = "ERROR: Node name " + new_name + " already exists\n";
            return false;
        }
        Node *node = node_access[old_name];
        node->setName(new_name);
        node_access.erase(old_name);
        node_access[new_name] = node;
        err = "SUCCESS: Node renamed from " + old_name + " to " + new_name + "\n";
        return true;
    }

    void BuildMNA() {
        nodetoindex.clear();
        ordernodes.clear();
        currentIndexmap.clear();
        int index = 0;
        for (Node *node: nodes) {
            if (!node->isGroundd()) {
                if (nodetoindex.count(node->getName()) == 0) {
                    nodetoindex[node->getName()] = index;
                    ordernodes.push_back(node);
                    index++;
                }
            }
        }
        int N = index;
        int M = 0;
        for (auto element: elements) {
            if (element->getCategory() == "Voltage Source" || element->getCategory() == "Controlled Voltage Source" || element->getCategory() == "Inductor")
                M++;
        }
        int size = M + N;
        vector<vector<double>> A(size, vector<double>(size, 0.0));
        vector<double> b(size, 0.0);
        int voltageindex = nodetoindex.size();
        for (auto element: elements) {
            if (element->getCategory() == "Voltage Source" || element->getCategory() == "Controlled Voltage Source" || element->getCategory() == "Inductor")
                currentIndexmap[element] = voltageindex++;
        }
        for (auto element: elements) {
            if (!element)
                continue;
            if(Analysistype=="DC"){
                if(auto *cap = dynamic_cast<Capacitor *>(element))
                    cap->DCstamp(A, b, nodetoindex, voltageindex, currentIndexmap);
                else if(auto *ind = dynamic_cast<Inductor *>(element))
                    ind->DCstamp(A, b, nodetoindex, voltageindex, currentIndexmap);
                else
                    element->stamp(A, b, nodetoindex, voltageindex, currentIndexmap);
            }
            else
                element->stamp(A, b, nodetoindex, voltageindex, currentIndexmap);
        }
        this->A = A;
        this->b = b;
    }

    void reset() {
        for (auto e: elements) delete e;
        for (auto n: nodes) delete n;

        elements.clear();
        nodes.clear();
        node_access.clear();
        nodetoindex.clear();
        ordernodes.clear();
        currentIndexmap.clear();
        gnd = nullptr;
    }

    vector<vector<double>> getMatrix() const { return A; }

    vector<double> getRHS() const { return b; }

    void clear () {
        for (auto* n : ordernodes) {
            n->setVoltage(0);
        }
        for (auto* e : elements) {
            e->setCurrent(0);
            if (auto* c = dynamic_cast<Capacitor*>(e)) c->updateVoltage(0);
            if (auto* l = dynamic_cast<Inductor*>(e)) l->updateCurrent(0);
        }
    }

};

class Resistor : public Element {
public:
    double voltage;

    Resistor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Resistor"; }


    double getVoltage() override {
        voltage = node1->getVoltage() - node2->getVoltage();
        return voltage;
    }

    double getCurrent() override {
        return (node1->getVoltage() - node2->getVoltage()) / value;
    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        int i = nodeToIndex[getFirstNode()->getName()];
        int j = nodeToIndex[getSecondNode()->getName()];
        double g = 1.0 / getValue();
        string n1 = node1->getName();
        string n2 = node2->getName();
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();

        if (hasN1 && hasN2) {
            int ip = nodeToIndex[n1];
            int jp = nodeToIndex[n2];
            A[ip][ip] += g;
            A[jp][jp] += g;
            A[ip][jp] -= g;
            A[jp][ip] -= g;
        } else if (hasN1) {
            int ip = nodeToIndex[n1];
            A[ip][ip] += g;
        } else if (hasN2) {
            int jp = nodeToIndex[n2];
            A[jp][jp] += g;
        }
    }
};

class Diode : public Element {
protected:
    string model;
public:
    Diode(string name, Node *n1, Node *n2) :
            Element(name, 0, n1, n2), model("D") {}

    string getType() override { return "Diode"; }

    virtual double calculateCurrent(double voltage) {

    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        int i = nodeToIndex[getFirstNode()->getName()];
        int j = nodeToIndex[getSecondNode()->getName()];
        double g = 1e9;

        A[i][i] += g;
        A[j][j] += g;
        A[i][j] -= g;
        A[j][i] -= g;
    }
};

class ZenerDiode : public Diode {
private:
    double breakdownvoltage;
public:
    ZenerDiode(string name, Node *n1, Node *n2, double breakdownvoltage_ = 5)
            : Diode(name, n1, n2), breakdownvoltage(breakdownvoltage_) {
        model = "Z";
    }

    string getType() override { return "Zener Diode"; }

    double calculateCurrent(double voltage) override {

    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
    }
};

class VoltageSource : public Element {
protected:
    double currenttime = 0;
public:
    VoltageSource(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    virtual double getVoltage(double time) = 0;

    string getType() override { return "Voltage Source"; }

    string getCategory() override { return "Voltage Source"; }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();
        int vIdx = currentIndexmap[this];
        if (hasN1) {
            int i = nodeToIndex[node1->getName()];
            A[i][vIdx] -= 1;
            A[vIdx][i] -= 1;
        }
        if (hasN2) {
            int j = nodeToIndex[node2->getName()];
            A[j][vIdx] += 1;
            A[vIdx][j] += 1;
        }
        b[vIdx] += getValue();
    }

    double getCurrent() override {
        return current;
    }
};

class DCVoltageSource : public VoltageSource {
public:
    DCVoltageSource(string name_, double value_, Node *n1, Node *n2)
            : VoltageSource(name_, value_, n1, n2) {}

    double getVoltage(double time) override {
        return getValue();
    }

    string getType() override { return "VDC"; }

    double getCurrent() override {
        return current;
    }
};

class SinusoidalVoltageSource : public VoltageSource {
private:
    double frequency;
    double offset;
    double amplitude;

public:
    SinusoidalVoltageSource(string name_, double amplitude_, double frequency_, double offset_, Node *n1, Node *n2)
            : VoltageSource(name_, amplitude_, n1, n2), amplitude(amplitude_), frequency(frequency_), offset(offset_) {}

    double getVoltage(double time) override {
        return getValue() * sin(2 * M_PI * frequency * time) + offset;
    }

    double getValue() override {
        return value * sin(2 * M_PI * frequency * currenttime) + offset;
    }

    double getFrequency() { return frequency; }

    double getOffset() { return offset; }

    double getAmp() {return amplitude;}

    string getType() override { return "VSIN"; }

    double getCurrent() override {
        return current;
    }

    void updateTime(double t) override { currenttime = t; }
};

class PulseVoltageSource : public VoltageSource {
private:
    double V1;
    double V2;
    double TD;
    double TR;
    double TF;
    double TOn;
    double period;

public:
    PulseVoltageSource(string name_, double V1_, double V2_, double TD_, double TR_, double TF_, double TOn_,
                       double period_, Node *n1, Node *n2)
            : VoltageSource(name_, V1_, n1, n2), V1(V1_), V2(V2_), TD(TD_), TR(TR_), TF(TF_), TOn(TOn_),
              period(period_) {}


    double getV1() const { return V1; }

    double getV2() const { return V2; }

    double getTD() const { return TD; }

    double getTR() const { return TR; }

    double getTF() const { return TF; }

    double getTOn() const { return TOn; }

    double getPeriod() const { return period; }

    double getVoltage(double time) override {
        double t = fmod(time, period);

        if (t < TD) {
            return V1;
        } else if (t < TD + TR) {
            return V1 + (V2 - V1) * (t - TD) / TR;
        } else if (t < TD + TR + TOn) {
            return V2;
        } else if (t < TD + TR + TOn + TF) {
            return V2 - (V2 - V1) * (t - TD - TOn - TR) / TF;
        } else {
            return V1;
        }
    }

    double getValue() override {
        return getVoltage(currenttime);
    }

    void updateTime(double t) { currenttime = t; }

    double getCurrent() override {
        return current;
    }

    string getType() override { return "VPULSE"; }
};

class CurrentSource : public Element {
public:
    CurrentSource(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    virtual double getCurrent(double time) = 0;

    string getType() override { return "Current Source"; }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexMap) override {
        bool hasN1 = nodeToIndex.count(node1->getName());
        bool hasN2 = nodeToIndex.count(node2->getName());
        double I = getValue();
        if (hasN1)
            b[nodeToIndex[node1->getName()]] -= I;
        if (hasN2)
            b[nodeToIndex[node2->getName()]] += I;

    }
};

class DCCurrentSource : public CurrentSource {
public:
    DCCurrentSource(string name_, double value_, Node *n1, Node *n2)
            : CurrentSource(name_, value_, n1, n2) {}

    double getCurrent(double time) override {
        return getValue();
    }

    string getType() override { return "IDC"; }

    double getCurrent() override {
        return current;
    }
};

class SinusoidalCurrentSource : public CurrentSource {
private:
    double frequency;
    double offset;
    double amplitude;

public:
    SinusoidalCurrentSource(string name_, double amplitude_, double frequency_, double offset_, Node *n1, Node *n2)
            : CurrentSource(name_, amplitude_, n1, n2), amplitude(amplitude_), frequency(frequency_), offset(offset_) {}

    double getCurrent(double time) override {
        return getValue() * sin(2 * M_PI * frequency * time) + offset;
    }

    double getFrequency() { return frequency; }

    double getOffset() { return offset; }

    double getValue() override {
        return value * sin(2 * M_PI * frequency * currenttime) + offset;
    }

    double getAmp() {return amplitude;}

    void updateTime(double t) override { currenttime = t; }

    string getType() override { return "ISIN"; }
};

class PulseCurrentSource : public CurrentSource {
private:
    double I1;
    double I2;
    double TD;
    double TR;
    double TF;
    double TOn;
    double period;

public:
    PulseCurrentSource(string name_, double I1_, double I2_, double TD_, double TR_, double TF_, double TOn_,
                       double period_, Node *n1, Node *n2)
            : CurrentSource(name_, I1_, n1, n2), I1(I1_), I2(I2_), TD(TD_), TR(TR_), TF(TF_), TOn(TOn_),
              period(period_) {}

    double getI1() const { return I1; }

    double getI2() const { return I2; }

    double getTD() const { return TD; }

    double getTR() const { return TR; }

    double getTF() const { return TF; }

    double getTOn() const { return TOn; }

    double getPeriod() const { return period; }

    double getCurrent(double time) override {
        double t = fmod(time, period);

        if (t < TD) {
            return I1;
        } else if (t < TD + TR) {
            return I1 + (I2 - I1) * (t - TD) / TR;
        } else if (t < TD + TR + TOn) {
            return I2;
        } else if (t < TD + TR + TOn + TF) {
            return I2 - (I2 - I1) * (t - TD - TOn - TR) / TF;
        } else {
            return I1;
        }
    }

    double getValue() override {
        return getCurrent(currenttime);
    }

    string getType() override { return "IPULSE"; }
};

class VCVS : public VoltageSource {
private:
    double gain;
    Node *controlNode1;
    Node *controlNode2;
public:
    VCVS(string name_, double gain_, Node *n1, Node *n2, Node *controlNode1_, Node *controlNode2_) :
            VoltageSource(name_, 0, n1, n2), gain(gain_), controlNode1(controlNode1_), controlNode2(controlNode2_) {}

    Node *getControlNode1() { return controlNode1; }

    Node *getControlNode2() { return controlNode2; }

    double getGain() { return gain; }

    string getCategory() override { return "Controlled Voltage Source"; }

    double getVoltage(double time) override {
        double controlVoltage = controlNode1->getVoltage() - controlNode2->getVoltage();
        return gain * controlVoltage;
    }

    string getType() override {
        return "VCVS";
    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        if (!currentIndexmap.count(this)) return;
        int vIdx = currentIndexmap[this];
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();
        bool hasN3 = !controlNode1->isGroundd();
        bool hasN4 = !controlNode2->isGroundd();

        int i = hasN1 ? nodeToIndex[node1->getName()] : -1;
        int j = hasN2 ? nodeToIndex[node2->getName()] : -1;
        int m = hasN3 ? nodeToIndex[controlNode1->getName()] : -1;
        int n = hasN4 ? nodeToIndex[controlNode2->getName()] : -1;

        if (hasN1) {
            A[i][vIdx] += 1;
            A[vIdx][i] += 1;
        }
        if (hasN2) {
            A[j][vIdx] -= 1;
            A[vIdx][j] -= 1;
        }
        if (hasN3 && hasN4) {
            A[vIdx][m] -= gain;
            A[vIdx][n] += gain;
        }
    }

    double getCurrent() override {
        return current;
    }
};

class CCVS : public VoltageSource {
private:
    double gain;
    Element *controlElement;
public:
    CCVS(string name_, double gain_, Node *n1, Node *n2, Element *controlElement_)
            : VoltageSource(name_, 0, n1, n2), gain(gain_), controlElement(controlElement_) {}

    double getGain() { return gain; }

    double getVoltage(double time) override {
        double controlCurrent = controlElement->getCurrent();
        return gain * controlCurrent;
    }

    Element *getcontrolElement() { return controlElement; }

    string getType() override {
        return "CCVS";
    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexMap) override {
        int vIdx = currentIndexMap[this];
        currentIndexMap[this] = vIdx;
        Element *control = getcontrolElement();
        int ctrlIdx = currentIndexMap[control];
        double gain = getGain();
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();
        int i = hasN1 ? nodeToIndex[node1->getName()] : -1;
        int j = hasN2 ? nodeToIndex[node2->getName()] : -1;

        if (hasN1) {
            A[i][vIdx] += 1;
            A[vIdx][i] += 1;
        }
        if (hasN2) {
            A[j][vIdx] -= 1;
            A[vIdx][j] -= 1;
        }

        A[vIdx][ctrlIdx] -= gain;
    }

    double getCurrent() override {
        return current;
    }

};

class VCCS : public CurrentSource {
private:
    double gain;
    Node *controlNode1;
    Node *controlNode2;
public:
    VCCS(string name_, double gain_, Node *n1, Node *n2, Node *controlNode1_, Node *controlNode2_) :
            CurrentSource(name_, 0, n1, n2), gain(gain_), controlNode1(controlNode1_), controlNode2(controlNode2_) {}

    Node *getControlNode1() { return controlNode1; }

    Node *getControlNode2() { return controlNode2; }

    double getGain() { return gain; }

    double getCurrent(double time) override {
        double controlVoltage = controlNode1->getVoltage() - controlNode2->getVoltage();
        return gain * controlVoltage;
    }

    string getCategory() override { return "Controlled Current Source"; }

    string getType() override {
        return "VCCS";
    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexmap) override {
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();
        bool hasN3 = !controlNode1->isGroundd();
        bool hasN4 = !controlNode2->isGroundd();

        int i = hasN1 ? nodeToIndex[node1->getName()] : -1;
        int j = hasN2 ? nodeToIndex[node2->getName()] : -1;
        int m = hasN3 ? nodeToIndex[controlNode1->getName()] : -1;
        int n = hasN4 ? nodeToIndex[controlNode2->getName()] : -1;

        if (hasN1 && hasN3) A[i][m] += gain;
        if (hasN1 && hasN4) A[i][n] -= gain;
        if (hasN2 && hasN3) A[j][m] -= gain;
        if (hasN2 && hasN4) A[j][n] += gain;
    }

    double getCurrent() override {
        return current;
    }
};

class CCCS : public CurrentSource {
private:
    double gain;
    Element *controlElement;
public:
    CCCS(string name_, double gain_, Node *n1, Node *n2, Element *controlElement_)
            : CurrentSource(name_, 0, n1, n2), gain(gain_), controlElement(controlElement_) {}

    double getGain() { return gain; }

    double getCurrent(double time) override {
        double controlCurrent = controlElement->getCurrent();
        return gain * controlCurrent;
    }

    Element *getcontrolElement() { return controlElement; }

    string getType() override {
        return "CCCS";
    }

    void stamp(vector<vector<double>> &A,
               vector<double> &b,
               map<string, int> &nodeToIndex,
               int &voltageIndex,
               map<Element *, int> &currentIndexMap) override {
        Element *control = getcontrolElement();

        int ctrlIdx = currentIndexMap[control];
        double gain = getGain();
        bool hasN1 = !node1->isGroundd();
        bool hasN2 = !node2->isGroundd();

        int i = hasN1 ? nodeToIndex[node1->getName()] : -1;
        int j = hasN2 ? nodeToIndex[node2->getName()] : -1;

        if (hasN1) A[i][ctrlIdx] += gain;
        if (hasN2) A[j][ctrlIdx] -= gain;
    }

    double getCurrent() override {
        return current;
    }
};

class Controller {
public:
    string handleError(string name, Circuit *circuit, double value) {
        string type;
        switch (name[0]) {
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                if (name.size() == 0)
                    type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        if (value <= 0 && (type == "Resistor" || type == "Capacitor" || type == "Inductor"))
            return "Error: " + type + " cannot be zero or negative\n";
        for (auto element: circuit->getElements()) {
            if (!element)
                continue;
            if (element->getName() == name) {
                switch (name[0]) {
                    case 'R':
                        return "Resistor " + name + " already exists in the circuit\n";
                    case 'C':
                        return "Capacitor " + name + " already exists in the circuit\n";
                    case 'L':
                        return "Inductor " + name + " already exists in the circuit\n";
                    case 'D':
                        return "Diode" + name + " already exists in the circuit\n";
                    default:
                        return "Source already exists in the circuit\n";
                }
            }
        }
        return "";
    }

    string addNewElement(string node1, string node2, string name, double value, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        switch (name[0]) {
            case 'R':
                element = new Resistor(name, value, n1, n2);
                break;
            case 'L':
                element = new Inductor(name, value, n1, n2);
                break;
            case 'C':
                element = new Capacitor(name, value, n1, n2);
                break;
        }
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string removeElement(string name, Circuit *circuit) {
        auto &elements = circuit->getElements();
        for (auto it = elements.begin(); it != elements.end(); it++) {
            if ((*it)->getName() == name) {
                (*it)->getFirstNode()->removeConnectedElement(*it);
                (*it)->getSecondNode()->removeConnectedElement(*it);
                Element *deleting = *it;
                elements.erase(it);
                delete deleting;
                return name + " removed successfully!\n";
            }
        }
        switch (name[0]) {
            case 'R':
                return "Error: Cannot delete resistor; component not found\n";
            case 'L':
                return "Error: Cannot delete inductor; component not found\n";
            case 'C':
                return "Error: Cannot delete capacitor; component not found\n";
            case 'D':
                return "Error: Cannot delete diode; component not found\n";
            default:
                return "Element not found!\n";
        }
    }

    string addDiode(string node1, string node2, string name, string model, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        if (model == "Z") {
            ZenerDiode *zenerDiode = new ZenerDiode(name, n1, n2);
            circuit->addElement(zenerDiode);
            n1->addConnectedElement(zenerDiode);
            n2->addConnectedElement(zenerDiode);
            return "Zener Diode added successfully!\n";
        } else if (model == "D") {
            Diode *diode = new Diode(name, n1, n2);
            circuit->addElement(diode);
            n1->addConnectedElement(diode);
            n2->addConnectedElement(diode);
            return "Diode added successfully!\n";
        }
    }

    string addGround(string node, Circuit *circuit) {
        if (circuit->getGround())
            return "Error! ground node already exists in the circuit\n";
        Node *n1 = circuit->getCreateNode(node);
        n1->setGround();
        circuit->setGround(node);
        return "Ground added to node " + node + " successfully!\n";

    }

    string deleteGround(string node, Circuit *circuit) {
        if (!circuit->getNodeAccess().count(node))
            return "Node does not exist\n";
        if (!circuit->removeGround(node))
            return "Error! Node is not set as ground\n";
        return "Ground removed from node " + node + " successfully!\n";
    }

    void showNodes(Circuit *circuit) {
        cout << "Available nodes:\n";
        if (!circuit->getNodeAccess().size())
            cout << "No node exists in the circuit\n";
        int index = 1;
        for (auto node: circuit->getNodeAccess()) {
            cout << node.first;
            if (index != circuit->getNodeAccess().size())
                cout << ", ";
            else
                cout << "\n";
            index++;
        }
    }

    void list(Circuit *circuit) {
        cout << "Elements:\n";
        if (!circuit->getElements().size())
            cout << "No elements exists in the circuit\n";
        for (auto element: circuit->getElements()) {
            cout << element->getType() << ": " << element->getName() << ", value: " << element->getValue() << endl;
        }
    }

    void listElement(string type, Circuit *circuit) {
        switch (type[0]) {
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        cout << type << "(s):\n";
        int index = 0;
        for (auto element: circuit->getElements()) {
            if (element->getType() == type) {
                cout << element->getType() << ": " << element->getName() << ", value: " << element->getValue() << endl;
                index++;
            }
        }
        if (!index)
            cout << "No " << type << " exist in the circuit\n";
    }

    string addDCVoltageSource(string name, string node1, string node2, double value, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new DCVoltageSource(name, value, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addSinusoidalVoltageSource(string name, string node1, string node2, double amplitude,
                                      double frequency, double offset, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new SinusoidalVoltageSource(name, amplitude, frequency, offset, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addPulseVoltageSource(string name, string node1, string node2, double V1, double V2, double TD, double TR,
                                 double TF, double TOn, double period, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new PulseVoltageSource(name, V1, V2, TD, TR, TF, TOn, period, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addDCCurrentSource(string name, string node1, string node2, double value, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new DCCurrentSource(name, value, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addSinusoidalCurrentSource(string name, string node1, string node2, double amplitude,
                                      double frequency, double offset, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new SinusoidalCurrentSource(name, amplitude, frequency, offset, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addPulseCurrentSource(string name, string node1, string node2, double I1, double I2, double TD, double TR,
                                 double TF, double TOn, double period, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        element = new PulseCurrentSource(name, I1, I2, TD, TR, TF, TOn, period, n1, n2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addVCVS(string name, string node1, string node2, string controlNode1, string controlNode2, double gain,
                   Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Node *cn1 = circuit->getCreateNode(controlNode1);
        Node *cn2 = circuit->getCreateNode(controlNode2);
        Element *element = nullptr;
        element = new VCVS(name, gain, n1, n2, cn1, cn2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addCCVS(string name, string node1, string node2, string controlElement, double gain, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        auto elements = circuit->getElements();
        Element *element = nullptr;
        Element *cElement = nullptr;
        for (auto &e: elements) {
            if (e->getName() == controlElement) {
                cElement = e;
                break;
            }
        }
        element = new CCVS(name, gain, n1, n2, cElement);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addVCCS(string name, string node1, string node2, string controlNode1, string controlNode2, double gain,
                   Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Node *cn1 = circuit->getCreateNode(controlNode1);
        Node *cn2 = circuit->getCreateNode(controlNode2);
        Element *element = nullptr;
        element = new VCCS(name, gain, n1, n2, cn1, cn2);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string addCCCS(string name, string node1, string node2, string controlElement, double gain, Circuit *circuit) {
        if (node1 == node2)
            return "The nodes can NOT be the same!\n";
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        auto elements = circuit->getElements();
        Element *element = nullptr;
        Element *cElement = nullptr;
        for (auto &e: elements) {
            if (e->getName() == controlElement) {
                cElement = e;
                break;
            }
        }
        element = new CCCS(name, gain, n1, n2, cElement);
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    void showCircuitDetails(Circuit *circuit) {
        cout << "================= Circuit Details =================\n";
        cout << "\nElements:\n";

        for (auto element: circuit->getElements()) {
            cout << "Type: " << element->getType()
                 << "  |  Name: " << element->getName()
                 << "  |  Nodes: (" << element->getFirstNode()->getName()
                 << ", " << element->getSecondNode()->getName() << ")";

            if (auto *vsin = dynamic_cast<SinusoidalVoltageSource *>(element)) {
                cout << "  |  Amplitude: " << vsin->getAmp()
                     << "  |  Frequency: " << vsin->getFrequency()
                     << "  |  Offset: " << vsin->getOffset();
            }

            else if (auto *vpulse = dynamic_cast<PulseVoltageSource *>(element)) {
                cout << "  |  V1: " << vpulse->getV1()
                     << "  |  V2: " << vpulse->getV2()
                     << "  |  TD: " << vpulse->getTD()
                     << "  |  TR: " << vpulse->getTR()
                     << "  |  TOn: " << vpulse->getTOn()
                     << "  |  TF: " << vpulse->getTF()
                     << "  |  Period: " << vpulse->getPeriod();
            }

            else if (auto *isin = dynamic_cast<SinusoidalCurrentSource *>(element)) {
                cout << "  |  Amplitude: " << isin->getAmp()
                     << "  |  Frequency: " << isin->getFrequency()
                     << "  |  Offset: " << isin->getOffset();
            }

            else if (auto *ipulse = dynamic_cast<PulseCurrentSource *>(element)) {
                cout << "  |  I1: " << ipulse->getI1()
                     << "  |  I2: " << ipulse->getI2()
                     << "  |  TD: " << ipulse->getTD()
                     << "  |  TR: " << ipulse->getTR()
                     << "  |  TOn: " << ipulse->getTOn()
                     << "  |  TF: " << ipulse->getTF()
                     << "  |  Period: " << ipulse->getPeriod();
            }

            else if (auto *vcvs = dynamic_cast<VCVS *>(element)) {
                cout << "  |  Gain: " << vcvs->getGain()
                     << "  |  Control Nodes: (" << vcvs->getControlNode1()->getName()
                     << ", " << vcvs->getControlNode2()->getName() << ")";
            }

            else if (auto *vccs = dynamic_cast<VCCS *>(element)) {
                cout << "  |  Gain: " << vccs->getGain()
                     << "  |  Control Nodes: (" << vccs->getControlNode1()->getName()
                     << ", " << vccs->getControlNode2()->getName() << ")";
            }

            else if (auto *cc = dynamic_cast<CCCS *>(element)) {
                cout << "  |  Gain: " << cc->getGain()
                     << "  |  Controlled Element: " << cc->getcontrolElement()->getName();
            } else if (auto *cc = dynamic_cast<CCVS *>(element)) {
                cout << "  |  Gain: " << cc->getGain()
                     << "  |  Controlled Element: " << cc->getcontrolElement()->getName();
            } else {
                cout << "  |  Value: " << element->getValue();
            }

            cout << endl;
        }

        cout << "Node Details:\n";
        for (auto node: circuit->getNodeAccess()) {
            cout << node.first << " : " << endl;
            if (node.second->getConnectedElements().size() == 0)
                cout << "No element is connected to this node!\n";
            else {
                for (auto element: node.second->getConnectedElements())
                    cout << element->getName() << endl;
            }
        }

        cout << "\nGND:\n";
        if (circuit->getGround())
            cout << circuit->getGround()->getName() << endl;
        else
            cout << "No ground is specified for this circuit\n";

        cout << "====================================================\n";
    }

    bool isConnected(Circuit *circuit) {
        if (circuit->getNodes().empty()) return true;

        unordered_map<string, vector<string>> adjacency;
        for (Element *e: circuit->getElements()) {
            if (!e) continue;
            string a = e->getFirstNode()->getName();
            string b = e->getSecondNode()->getName();
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }

        unordered_set<string> visited;
        queue<string> q;
        string start = circuit->getNodes()[0]->getName();
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            string current = q.front();
            q.pop();
            for (const string &neighbor: adjacency[current]) {
                if (!visited.count(neighbor)) {
                    visited.insert(neighbor);
                    q.push(neighbor);
                }
            }
        }

        return visited.size() == circuit->getNodeAccess().size();
    }

    bool preAnalysisErrs(Circuit *circuit) {
        if (!circuit->getGround()) {
            cout << "Error: No ground node detected in the circuit.\n";
            return false;
        }
        if (!isConnected(circuit)) {
            cout << "Error: Not all nodes are connected!\n";
            return false;
        }
        for (auto element: circuit->getElements()) {
            if (element->getType() == "CCCS") {
                auto *cc = dynamic_cast<CCCS *>(element);
                if (cc) {
                    Element *ctrl = cc->getcontrolElement();
                    bool found = false;
                    for (auto e: circuit->getElements()) {
                        if (e == ctrl) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        cout << "Error: Dependent source " << cc->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }
            if (element->getType() == "CCVS") {
                auto *cv = dynamic_cast<CCVS *>(element);
                if (cv) {
                    Element *ctrl = cv->getcontrolElement();
                    bool found = false;
                    for (auto e: circuit->getElements()) {
                        if (e == ctrl) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        cout << "Error: Dependent source " << cv->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }
            if (element->getType() == "VCCS") {
                auto *vc = dynamic_cast<VCCS *>(element);
                if (vc) {
                    Node *ctrl = vc->getControlNode1();
                    Node *ctrl2 = vc->getControlNode2();
                    bool found = false;
                    bool found2 = false;
                    for (auto e: circuit->getNodes()) {
                        if (e == ctrl) {
                            found = true;
                        }
                        if (e == ctrl2)
                            found2 = true;
                        if (found && found2)
                            break;
                    }
                    if (!found || !found2 || (ctrl->getConnectedElements().size() == 0) ||
                        (ctrl2->getConnectedElements().size() == 0)) {
                        cout << "Error: Dependent source " << vc->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }
            if (element->getType() == "VCVS") {
                auto *vv = dynamic_cast<VCVS *>(element);
                if (vv) {
                    Node *ctrl = vv->getControlNode1();
                    Node *ctrl2 = vv->getControlNode2();
                    bool found = false;
                    bool found2 = false;
                    for (auto e: circuit->getNodes()) {
                        if (e == ctrl) {
                            found = true;
                        }
                        if (e == ctrl2)
                            found2 = true;
                        if (found && found2)
                            break;
                    }
                    if (!found || !found2 || (ctrl->getConnectedElements().size() == 0) ||
                        (ctrl2->getConnectedElements().size() == 0)) {
                        cout << "Error: Dependent source " << vv->getName()
                             << "  has an undefined control element.\n";
                        return false;
                    }
                }
            }
        }
        return true;
    }

    void loadCircuit(Circuit *circuit, string address) {
        ifstream fin(address, ios::in);

        if (!fin) {
            cout << "Error opening file!" << endl;
            return;
        }

        string line;
        regex word_regex("\\S+");
        string tmpAddress;
        unsigned int t1;
        unsigned int t2;
        t1 = address.find_last_of("\\") + 1;
        t2 = address.length() - 11;
        tmpAddress = address.substr(t1, t2);
        cout << tmpAddress << ":" << endl;
        while (getline(fin, line)) {
            auto words_begin = sregex_iterator(line.begin(), line.end(), word_regex);
            auto words_end = sregex_iterator();

            vector<string> words;

            for (auto it = words_begin; it != words_end; ++it) {
                words.push_back(it->str());
            }

            string name;
            string node1;
            string node2;
            string tmpValue;
            string tmpAmplitude;
            string tmpFrequency;
            string tmpOffset;
            string tmpVi;
            string tmpVf;
            string tmpTD;
            string tmpTR;
            string tmpTF;
            string tmpTOn;
            string tmpPeriod;


            string type = words[0];
            if (type == "R" || type == "L" || type == "C" || type == "I" || type == "V") {
                name = words[1];
                node1 = words[2];
                node2 = words[3];
                tmpValue = words[4];
                if (node1 == "GND") {
                    addGround(node1, circuit);
                }
                if (node2 == "GND") {
                    addGround(node2, circuit);
                }

                double value = 0;

                if (tmpValue[0] == '-') {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpValue.back() == 'G') {
                        value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        value = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (type == "R" || type == "C" || type == "L") {
                    cout << addNewElement(node1, node2, name, value, circuit);

                } else if (type == "V") {
                    cout << addDCVoltageSource(name, node1, node2, value, circuit);
                } else if (type == "I") {
                    cout << addDCCurrentSource(name, node1, node2, value, circuit);
                }
            }
            if (type == "D" || type == "Z") {
                name = words[1];
                node1 = words[2];
                node2 = words[3];
                if (node1 == "GND") {
                    addGround(node1, circuit);
                }
                if (node2 == "GND") {
                    addGround(node2, circuit);
                }
                string model;
                if (type == "D") {
                    model = "D";
                }
                if (type == "Z") {
                    model = "Z";
                }
                cout << addDiode(node1, node2, name, model, circuit);
            }
            if (type == "VSIN" || type == "ISIN") {
                name = words[1];
                node1 = words[2];
                node2 = words[3];
                tmpAmplitude = words[4];
                tmpFrequency = words[5];
                tmpOffset = words[6];
                double amplitude;
                double frequency;
                double offset;
                try {
                    if (tmpAmplitude.back() == 'G') {
                        amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e9;
                    } else if (tmpAmplitude.back() == 'M') {
                        amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e6;
                    } else if (tmpAmplitude.back() == 'k' || tmpAmplitude.back() == 'K') {
                        amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e3;
                    } else if (tmpAmplitude.back() == 'u') {
                        amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-6;
                    } else if (tmpAmplitude.back() == 'n') {
                        amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-9;
                    } else if (tmpAmplitude.back() == 'm') {
                        amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-3;
                    } else {
                        amplitude = stod(tmpAmplitude);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpFrequency.back() == 'G') {
                        frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e9;
                    } else if (tmpFrequency.back() == 'M') {
                        frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e6;
                    } else if (tmpFrequency.back() == 'k' || tmpFrequency.back() == 'K') {
                        frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e3;
                    } else if (tmpFrequency.back() == 'u') {
                        frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-6;
                    } else if (tmpFrequency.back() == 'n') {
                        frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-9;
                    } else if (tmpFrequency.back() == 'm') {
                        frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-3;
                    } else {
                        frequency = stod(tmpFrequency);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                try {
                    if (tmpOffset.back() == 'G') {
                        offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e9;
                    } else if (tmpOffset.back() == 'M') {
                        offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e6;
                    } else if (tmpOffset.back() == 'k' || tmpOffset.back() == 'K') {
                        offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e3;
                    } else if (tmpOffset.back() == 'u') {
                        offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-6;
                    } else if (tmpOffset.back() == 'n') {
                        offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-9;
                    } else if (tmpOffset.back() == 'm') {
                        offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-3;
                    } else {
                        offset = stod(tmpOffset);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                if (type == "VSIN") {
                    cout << addSinusoidalVoltageSource(name, node1, node2, amplitude, frequency,
                                                       offset, circuit);
                }
                if (type == "ISIN") {
                    cout << addSinusoidalCurrentSource(name, node1, node2, amplitude, frequency,
                                                       offset, circuit);
                }
            }

            if (type == "VPULSE" || type == "IPULSE") {
                name = words[1];
                node1 = words[2];
                node2 = words[3];
                tmpVi = words[4];
                tmpVf = words[5];
                tmpTD = words[6];
                tmpTR = words[7];
                tmpTF = words[8];
                tmpTOn = words[9];
                tmpPeriod = words[10];
                double Vi, Vf, TD, TR, TF, TOn, period;
                try {
                    if (tmpVi.back() == 'G') {
                        Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e9;
                    } else if (tmpVi.back() == 'M') {
                        Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e6;
                    } else if (tmpVi.back() == 'k' || tmpVi.back() == 'K') {
                        Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e3;
                    } else if (tmpVi.back() == 'u') {
                        Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-6;
                    } else if (tmpVi.back() == 'n') {
                        Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-9;
                    } else if (tmpVi.back() == 'm') {
                        Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-3;
                    } else {
                        Vi = stod(tmpVi);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpVf.back() == 'G') {
                        Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e9;
                    } else if (tmpVf.back() == 'M') {
                        Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e6;
                    } else if (tmpVf.back() == 'k' || tmpVf.back() == 'K') {
                        Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e3;
                    } else if (tmpVf.back() == 'u') {
                        Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-6;
                    } else if (tmpVf.back() == 'n') {
                        Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-9;
                    } else if (tmpVf.back() == 'm') {
                        Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-3;
                    } else {
                        Vf = stod(tmpVf);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpTD.back() == 'G') {
                        TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e9;
                    } else if (tmpTD.back() == 'M') {
                        TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e6;
                    } else if (tmpTD.back() == 'k' || tmpTD.back() == 'K') {
                        TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e3;
                    } else if (tmpTD.back() == 'u') {
                        TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-6;
                    } else if (tmpTD.back() == 'n') {
                        TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-9;
                    } else if (tmpTD.back() == 'm') {
                        TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-3;
                    } else {
                        TD = stod(tmpTD);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpTR.back() == 'G') {
                        TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e9;
                    } else if (tmpTR.back() == 'M') {
                        TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e6;
                    } else if (tmpTR.back() == 'k' || tmpTR.back() == 'K') {
                        TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e3;
                    } else if (tmpTR.back() == 'u') {
                        TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-6;
                    } else if (tmpTR.back() == 'n') {
                        TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-9;
                    } else if (tmpTR.back() == 'm') {
                        TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-3;
                    } else {
                        TR = stod(tmpTR);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpTF.back() == 'G') {
                        TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e9;
                    } else if (tmpTF.back() == 'M') {
                        TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e6;
                    } else if (tmpTF.back() == 'k' || tmpTF.back() == 'K') {
                        TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e3;
                    } else if (tmpTF.back() == 'u') {
                        TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-6;
                    } else if (tmpTF.back() == 'n') {
                        TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-9;
                    } else if (tmpTF.back() == 'm') {
                        TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-3;
                    } else {
                        TF = stod(tmpTF);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpTOn.back() == 'G') {
                        TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e9;
                    } else if (tmpTOn.back() == 'M') {
                        TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e6;
                    } else if (tmpTOn.back() == 'k' || tmpTOn.back() == 'K') {
                        TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e3;
                    } else if (tmpTOn.back() == 'u') {
                        TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-6;
                    } else if (tmpTOn.back() == 'n') {
                        TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-9;
                    } else if (tmpTOn.back() == 'm') {
                        TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-3;
                    } else {
                        TOn = stod(tmpTOn);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }

                try {
                    if (tmpPeriod.back() == 'G') {
                        period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e9;
                    } else if (tmpPeriod.back() == 'M') {
                        period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e6;
                    } else if (tmpPeriod.back() == 'k' || tmpPeriod.back() == 'K') {
                        period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e3;
                    } else if (tmpPeriod.back() == 'u') {
                        period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-6;
                    } else if (tmpPeriod.back() == 'n') {
                        period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-9;
                    } else if (tmpPeriod.back() == 'm') {
                        period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-3;
                    } else {
                        period = stod(tmpPeriod);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (type == "VPULSE") {
                    cout << addPulseVoltageSource(name, node1, node2, Vi, Vf, TD, TR, TF, TOn,
                                                  period, circuit);
                }
                if (type == "IPULSE") {
                    cout << addPulseCurrentSource(name, node1, node2, Vi, Vf, TD, TR, TF, TOn,
                                                  period, circuit);
                }
            }
            if (type == "E" || type == "G") {
                name = words[1];
                node1 = words[2];
                node2 = words[3];
                string controlNode1 = words[4];
                string controlNode2 = words[5];
                string tmpGain = words[5];
                double gain;
                Node *cn1 = circuit->getCreateNode(controlNode1);
                Node *cn2 = circuit->getCreateNode(controlNode2);
                bool cn1Connected = false;
                bool cn2Connected = false;
                for (auto element: circuit->getElements()) {
                    if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                        cn1Connected = true;
                    }
                    if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                        cn2Connected = true;
                    }
                }

                if (!cn1Connected || !cn2Connected) {
                    cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected"
                         << endl;
                    continue;
                }
                try {
                    if (tmpGain.back() == 'G') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e9;
                    } else if (tmpGain.back() == 'M') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e6;
                    } else if (tmpGain.back() == 'k' || tmpGain.back() == 'K') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e3;
                    } else if (tmpGain.back() == 'u') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-6;
                    } else if (tmpGain.back() == 'n') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-9;
                    } else if (tmpGain.back() == 'm') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpGain);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }
                if (type == "E") {
                    cout << addVCVS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                }
                if (type == "G") {
                    cout << addVCCS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                }
            }

            if (type == "H" || type == "F") {
                name = words[1];
                node1 = words[2];
                node2 = words[3];
                string cElement = words[4];
                string tmpGain = words[5];
                double gain;

                Element *controlElement = nullptr;
                for (auto &element: circuit->getElements()) {
                    if (element->getName() == cElement) {
                        controlElement = element;
                        break;
                    }
                }
                if (controlElement == nullptr) {
                    cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                    continue;
                }
                try {
                    if (tmpGain.back() == 'G') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e9;
                    } else if (tmpGain.back() == 'M') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e6;
                    } else if (tmpGain.back() == 'k' || tmpGain.back() == 'K') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e3;
                    } else if (tmpGain.back() == 'u') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-6;
                    } else if (tmpGain.back() == 'n') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-9;
                    } else if (tmpGain.back() == 'm') {
                        gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpGain);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }
                if (type == "H") {
                    cout << addCCVS(name, node1, node2, cElement, gain, circuit);
                }
                if (type == "F") {
                    cout << addCCCS(name, node1, node2, cElement, gain, circuit);
                }

            }

        }

        fin.close();
    }

    void saveCircuitToFile(Circuit *circuit, string filename) {
        const char *dir = "drafts";
        struct stat sb;
        if (stat(dir, &sb) == -1) {
            if (_mkdir(dir) == -1) {
                cout << "Error: Unable to create the 'drafts' folder.\n";
            }
        }

        if (filename.find(".txt") == string::npos)
            filename += ".txt";
        string filePath = "drafts\\" + filename;

        ofstream fout(filePath);
        if (!fout.is_open()) {
            cout << "Error: Could not open file '" + filePath + "' for writing.\n";
        }

        for (Element *e: circuit->getElements()) {
            if (!e) continue;
            string type = e->getType();
            string name = e->getName();
            string n1 = e->getFirstNode()->getName();
            string n2 = e->getSecondNode()->getName();

            if (type == "Resistor" || type == "Capacitor" || type == "Inductor") {
                double value = e->getValue();
                char prefix = (type == "Resistor" ? 'R' : type == "Capacitor" ? 'C' : 'L');
                fout << prefix << " " << name << " " << n1 << " " << n2 << " " << value << "\n";
            } else if (type == "VDC") {
                auto *vdc = dynamic_cast<DCVoltageSource *>(e);
                if (vdc) {
                    double value = vdc->getValue();
                    fout << "V " << name << " " << n1 << " " << n2 << " " << value << "\n";
                }
            } else if (type == "IDC") {
                auto *idc = dynamic_cast<DCCurrentSource *>(e);
                if (idc) {
                    double value = idc->getValue();
                    fout << "I " << name << " " << n1 << " " << n2 << " " << value << "\n";
                }
            } else if (type == "Diode") {
                fout << "D " << name << " " << n1 << " " << n2 << "\n";
            } else if (type == "Zener Diode") {
                fout << "Z " << name << " " << n1 << " " << n2 << "\n";
            } else if (type == "VSIN") {
                auto *sv = dynamic_cast<SinusoidalVoltageSource *>(e);
                if (sv) {
                    double amp = sv->getValue();
                    double freq = sv->getFrequency();
                    double off = sv->getOffset();
                    fout << "VSIN " << name << " " << n1 << " " << n2
                         << " " << amp << " " << freq << " " << off << "\n";
                }
            } else if (type == "ISIN") {
                auto *sc = dynamic_cast<SinusoidalCurrentSource *>(e);
                if (sc) {
                    double amp = sc->getValue();
                    double freq = sc->getFrequency();
                    double off = sc->getOffset();
                    fout << "ISIN " << name << " " << n1 << " " << n2
                         << " " << amp << " " << freq << " " << off << "\n";
                }
            } else if (type == "VPULSE") {
                auto *pv = dynamic_cast<PulseVoltageSource *>(e);
                if (pv) {
                    double V1 = pv->getV1();
                    double V2 = pv->getV2();
                    double TD = pv->getTD();
                    double TR = pv->getTR();
                    double TF = pv->getTF();
                    double TOn = pv->getTOn();
                    double period = pv->getPeriod();
                    fout << "VPULSE " << name << " " << n1 << " " << n2
                         << " " << V1 << " " << V2 << " " << TD << " " << TR
                         << " " << TF << " " << TOn << " " << period << "\n";
                }
            } else if (type == "IPULSE") {
                auto *pc = dynamic_cast<PulseCurrentSource *>(e);
                if (pc) {
                    double I1 = pc->getI1();
                    double I2 = pc->getI2();
                    double TD = pc->getTD();
                    double TR = pc->getTR();
                    double TF = pc->getTF();
                    double TOn = pc->getTOn();
                    double period = pc->getPeriod();
                    fout << "IPULSE " << name << " " << n1 << " " << n2
                         << " " << I1 << " " << I2 << " " << TD << " " << TR
                         << " " << TF << " " << TOn << " " << period << "\n";
                }
            } else if (type == "VCVS") {
                auto *vd = dynamic_cast<VCVS *>(e);
                if (vd) {
                    string cn1 = vd->getControlNode1()->getName();
                    string cn2 = vd->getControlNode2()->getName();
                    double gain = vd->getGain();
                    fout << "E " << name << " " << n1 << " " << n2
                         << " " << cn1 << " " << cn2 << " " << gain << "\n";
                }
            } else if (type == "VCCS") {
                auto *vc = dynamic_cast<VCCS *>(e);
                if (vc) {
                    string cn1 = vc->getControlNode1()->getName();
                    string cn2 = vc->getControlNode2()->getName();
                    double gain = vc->getGain();
                    fout << "G " << name << " " << n1 << " " << n2
                         << " " << cn1 << " " << cn2 << " " << gain << "\n";
                }
            } else if (type == "CCVS") {
                auto *cv = dynamic_cast<CCVS *>(e);
                if (cv) {
                    string ctrlName = cv->getcontrolElement()->getName();
                    double gain = cv->getGain();
                    fout << "H " << name << " " << n1 << " " << n2
                         << " " << ctrlName << " " << gain << "\n";
                }
            } else if (type == "CCCS") {
                auto *cc = dynamic_cast<CCCS *>(e);
                if (cc) {
                    string ctrlName = cc->getcontrolElement()->getName();
                    double gain = cc->getGain();
                    fout << "F " << name << " " << n1 << " " << n2
                         << " " << ctrlName << " " << gain << "\n";
                }
            }
        }

        fout.close();
    }

    void showAndLoadSchematic(Circuit *circuit, string &loadedFilename) {
        const char *dirname = "drafts";
        DIR *dirp = opendir(dirname);
        if (!dirp) {
            cout << "Error: Could not open the 'drafts' folder\n";
        } else {
            vector<string> files;
            struct dirent *entry;
            while ((entry = readdir(dirp)) != nullptr) {
                string name = entry->d_name;
                if (name == "." || name == "..")
                    continue;
                files.push_back(name);
            }
            closedir(dirp);

            if (files.empty()) {
                cout << "No existing schematics found in 'drafts'.\n";
            } else {
                while (true) {
                    cout << "-choose existing schematic:\n";
                    for (int i = 0; i < files.size(); ++i) {
                        cout << "  " << (i + 1) << "-" << files[i] << "\n";
                    }

                    cout << "Type 'return' to exit the menu.\n";

                    int choice;
                    string user_input;

                    cout << "Enter your choice: ";
                    getline(cin, user_input);

                    if (user_input == "return") {
                        cout << "Continue building your circuit" << endl;
                        break;
                    }

                    try {
                        choice = stoi(user_input);
                        if (choice < 1 || choice > (int) files.size()) {
                            cout << "-Error: Inappropriate input\n";
                        } else {
                            string selectedName = files[choice - 1];
                            string fullPath = string(dirname) + "\\" + selectedName;
                            loadCircuit(circuit, fullPath);
                            loadedFilename = selectedName;
                            break;
                        }
                    } catch (const invalid_argument& e) {
                        cout << "-Error: Inappropriate input\n";
                    }
                }
            }
        }
    }

    void checkMatrix(Circuit *circuit) {
        circuit->BuildMNA();
        const auto &A = circuit->getMatrix();
        const auto &b = circuit->getRHS();
        int n = A.size();
        if (n == 0) {
            cout << "MNA matrix is empty!\n";
            return;
        }

        cout << "MNA Matrix A:\n";
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < A[i].size(); ++j) {
                cout << setw(10) << fixed << setprecision(3) << A[i][j] << " ";
            }
            cout << "\n";
        }

        cout << "\nVector b:\n";
        for (int i = 0; i < b.size(); ++i) {
            cout << "b[" << i << "] = " << fixed << setprecision(3) << b[i] << "\n";
        }
        cout << "\n";
    }

    void DCAnalysis(Circuit *circuit) {
        circuit->setType("DC");
        circuit->clear();
        circuit->BuildMNA();
        vector<vector<double>> A = circuit->getMatrix();
        vector<double> b = circuit->getRHS();
        vector<double> solution = gaussianElimination(A, b);
        for (int i = 0; i < circuit->ordernodes.size(); ++i) {
            circuit->ordernodes[i]->setVoltage(solution[i]);
        }
        for (Element *e: circuit->getElements()) {
            if (circuit->currentIndexmap.count(e)) {
                int idx = circuit->currentIndexmap[e];
                e->setCurrent(solution[idx]);
            }
        }
        cout << "========== DC Operating Point Analysis ==========\n";
        cout << "\nNode Voltages:\n";
        for (auto *node: circuit->ordernodes) {
            cout << "V(" << node->getName() << ") = " << fixed << setprecision(3) << node->getVoltage() << " V\n";
        }
        cout << "\nElement Currents:\n";
        for (Element *e: circuit->getElements()) {
            cout << "I(" << e->getName() << ") = " << fixed << setprecision(3) << e->getCurrent() << " A\n";
        }

        cout << "=================================================\n";
    }

    void DCSweep(string sweptSource, double startValue, double endValue, double increment,
                 string type, string targetName, Circuit *circuit) {
        circuit->setType("DC");
        circuit->clear();
        Element *sweepElem = nullptr;
        for (auto *e: circuit->getElements()) {
            if (e->getName() == sweptSource) {
                sweepElem = e;
                break;
            }
        }

        if (!sweepElem) {
            cout << "Error: source " << sweptSource << " not found!\n";
            return;
        }

        if (type == "V" && !circuit->getNodeAccess().count(targetName)) {
            cout << "Error: Node " << targetName << " not found in the circuit\n";
            return;
        } else if (type == "I") {
            bool found = false;
            for (auto &element: circuit->getElements()) {
                if (element->getName() == targetName) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout << "Error: Element " << targetName << " not found in the circuit\n";
                return;
            }
        }

        double originalValue = sweepElem->getValue();
        cout << "=============== DC Sweep Analysis ===============\n";

        for (double val = startValue; val <= endValue + 1e-9; val += increment) {
            sweepElem->setValue(val);
            circuit->BuildMNA();
            vector<vector<double>> A = circuit->getMatrix();
            vector<double> b = circuit->getRHS();
            vector<double> solution = gaussianElimination(A, b);

            for (int i = 0; i < circuit->ordernodes.size(); ++i)
                circuit->ordernodes[i]->setVoltage(solution[i]);

            int N = circuit->ordernodes.size();

            for (Element *e: circuit->getElements()) {
                if (circuit->currentIndexmap.count(e)) {
                    int idx = circuit->currentIndexmap[e];
                    double current = solution[idx];
                    e->setCurrent(current);
                }
            }

            cout << sweptSource << " = " << val << " : ";

            if (type == "V") {
                if (!circuit->getNodeAccess().count(targetName)) {
                    cout << "Node " << targetName << " not found in the circuit" << endl;
                    return;
                }
                bool found = false;
                for (auto *n: circuit->ordernodes) {
                    if (n->getName() == targetName) {
                        cout << "V(" << targetName << ") = " << fixed << setprecision(3) << n->getVoltage() << " V";
                        found = true;
                        break;
                    }
                }
                if (!found) cout << "Error: node " << targetName << " not found";
            } else if (type == "I") {
                Element *analyzedElement = nullptr;
                bool found = false;
                for (auto &element: circuit->getElements()) {
                    if (element->getName() == targetName) {
                        analyzedElement = element;
                        break;
                    }
                }

                if (!analyzedElement) {
                    cout << "Component " << targetName << " not found in the circuit" << endl;
                    return;
                }

                for (auto *e: circuit->getElements()) {
                    if (e->getName() == targetName) {
                        cout << "I(" << targetName << ") = " << fixed << setprecision(3) << e->getCurrent() << " A";
                        found = true;
                        break;
                    }
                }
                if (!found) cout << "Error: element " << targetName << " not found";
            }
            cout << endl;
        }
        cout << "=================================================\n";

        sweepElem->setValue(originalValue);
    }

    void multipleDCSweep(vector<string> &probes, string sweptSource, double startValue, double endValue, double increment, Circuit* circuit){
        for (auto& p : probes){
            if (p[0]=='V'){
                string type = "V";
                p.erase(0, 1);
                p.erase(remove(p.begin(), p.end(), '('), p.end());
                p.erase(remove(p.begin(), p.end(), ')'), p.end());
                DCSweep(sweptSource, startValue, endValue, increment, type, p, circuit);
            }
            if (p[0] == 'I'){
                string type = "I";
                p.erase(0, 1);
                p.erase(remove(p.begin(), p.end(), '('), p.end());
                p.erase(remove(p.begin(), p.end(), ')'), p.end());
                DCSweep(sweptSource, startValue, endValue, increment, type, p, circuit);
            }
        }
    }


    void TransientAnalysis(double startValue, double endValue, double increment,
                           string type, string target, Circuit *circuit) {
        circuit->setType("Transient");
        circuit->clear();
        if (type == "V" && !circuit->getNodeAccess().count(target)) {
            cout << "Error: Node " << target << " not found in the circuit\n";
            return;
        } else if (type == "I") {
            bool found = false;
            for (auto &element: circuit->getElements()) {
                if (element->getName() == target) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                cout << "Error: Element " << target << " not found in the circuit\n";
                return;
            }
        }

        double delta_t = 1e-3;
        circuit->setDeltat(delta_t);

        int totalSteps = (int)((endValue - startValue) / delta_t) + 1;
        double nextPrintTime = startValue;

        cout << "========== Transient Analysis ==========\n";
        cout << "Analyzing: " << (type == "V" ? "V(" + target + ")" : "I(" + target + ")") << endl;

        for (int step = 0; step <= totalSteps; ++step) {
            double t = startValue + step * delta_t;

            for (Element* e : circuit->getElements()) {
                e->updateTime(t);
            }
            circuit->BuildMNA();
            auto A = circuit->getMatrix();
            auto b = circuit->getRHS();
            auto solution = gaussianElimination(A, b);

            for (int i = 0; i < circuit->ordernodes.size(); ++i)
                circuit->ordernodes[i]->setVoltage(solution[i]);

            for (Element* e : circuit->getElements()) {
                if (circuit->currentIndexmap.count(e)) {
                    double inew=solution[circuit->currentIndexmap[e]];
                    e->setCurrent(inew);
                    if (auto* ind = dynamic_cast<Inductor*>(e)) {
                        ind->updateCurrent(inew);
                    }
                }
            }

            for (Element* e : circuit->getElements()) {
                if (auto* cap = dynamic_cast<Capacitor*>(e)) {
                    double vcap = cap->getFirstNode()->getVoltage() - cap->getSecondNode()->getVoltage();
                    cap->updateVoltage(vcap);
                }
            }



            if (abs(t - nextPrintTime) < delta_t / 2 || t >= endValue) {
                cout << fixed << setprecision(3);
                cout << "t = " << t << " s : ";

                if (type == "V") {
                    for (auto* n : circuit->ordernodes) {
                        if (n->getName() == target) {
                            cout << "V(" << target << ") = " << clean(n->getVoltage()) << " V";
                            break;
                        }
                    }
                } else if (type == "I") {
                    for (auto* e : circuit->getElements()) {
                        if (e->getName() == target) {
                            cout << "I(" << target << ") = " << clean(e->getCurrent()) << " A";
                            break;
                        }
                    }
                }

                cout << endl;
                nextPrintTime += increment;
            }
        }

        cout << "========================================\n";
    }

    void multipleTransient(vector<string> &probes, double startValue, double endValue, double increment, Circuit* circuit) {
        for (auto& p : probes) {
            if (p[0] == 'V') {
                string type = "V";
                p.erase(0, 1);
                p.erase(remove(p.begin(), p.end(), '('), p.end());
                p.erase(remove(p.begin(), p.end(), ')'), p.end());
                TransientAnalysis(startValue, endValue, increment, type, p, circuit);
            }
            if (p[0] == 'I') {
                string type = "I";
                p.erase(0, 1);
                p.erase(remove(p.begin(), p.end(), '('), p.end());
                p.erase(remove(p.begin(), p.end(), ')'), p.end());
                TransientAnalysis(startValue, endValue, increment, type, p, circuit);
            }
        }
    }



};


class View {
private:
    Circuit *circuit;
    Controller controller;
    string loadedFilename;
public:
    void run() {
        circuit = new Circuit();
        string input;
        smatch match;
        regex add_element(
                R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex remove_element(R"(^\s*delete\s+([A-Za-z])(\w+)\s*$)");
        regex add_diode(R"(^\s*add\s+(D)(\w+)\s+(\S+)\s+(\S+)\s+(D|Z)\s*$)");
        regex show_details(R"(^\s*show\s+details\s*$)");
        regex add_ground(R"(^\s*add\s+(\w+)\s+(\w+)\s*$)");
        regex delete_ground(R"(^\s*delete\s+(\w+)\s+(\w+)\s*$)");
        regex show_nodes(R"(^\s*.nodes\s*$)");
        regex list(R"(^\s*.list\s*$)");
        regex list_element(R"(^\s*.list\s+(\w+)\s*$)");
        regex rename_node(R"(^\s*\.rename\s+node\s+(\S+)\s+(\S+)\s*$)");
        regex new_file(R"(NewFile ([a-zA-Z0-9\\-_\\.:\\/()\\s]+))");
        regex addDCSource(
                R"(add\s+(VoltageSource|CurrentSource)(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex addSINSource(
                R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+SIN\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex addPULSESource(
                R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+PULSE\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex addVCVS(
                R"(^\s*add\s+E(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex addVCCS(
                R"(^\s*add\s+G(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex addCCVS(R"(^\s*add\s+H(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex addCCCS(R"(^\s*add\s+F(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?[GMKkmnp]?)\s*$)");
        regex preanalysischeck(R"(^\s*pre-analysis\s+check\s*$)");
        regex save_file(R"(^\s*save\s+file\s+(\w+)\s*$)");
        regex show_existing_schematics(R"(^\s*show\s+existing\s+schematics\s*$)");
        regex check_matrix(R"(^\s*check\s+matrix\s*$)");
        regex DC_Analaysis(R"(^\s*DC\s+op\s+analysis\s*$)");
        regex DCSweep(R"(^\s*\.DC\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([VI])\((\S+)\)\s*$)", regex::icase);
        regex multipleDCSweep(R"(^\.DC\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+[IV]\(\S+\)){2,}\s*$)", regex::icase);
        regex multipleTRAN(R"(^\.TRAN\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+[IV]\(\S+\)){2,}\s*$)", regex::icase);
        regex transient(R"(^\s*\.TRAN\s+(\S+)\s+(\S+)\s+(\S+)\s+([VI])\((\S+)\)\s*$)", regex::icase);
        regex reset(R"(^reset)");
        regex exit(R"(^exit$)");
        while (true) {
            getline(cin, input);
            if (regex_match(input, match, addDCSource)) {
                string name = match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number = match[5].str();
                string unit = match[6].str();
                double value;
                try {
                    value = stod(number);
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                string error = controller.handleError(name, circuit, value);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (!unit.empty()) {
                    switch (unit[0]) {
                        case 'G':
                            value *= 1e9;
                            break;
                        case 'M':
                            value *= 1e6;
                            break;
                        case 'k':
                        case 'K':
                            value *= 1e3;
                            break;
                        case 'm':
                            value *= 1e-3;
                            break;
                        case 'u':
                            value *= 1e-6;
                            break;
                        case 'n':
                            value *= 1e-9;
                            break;
                    }
                }
                if (match[1] == "VoltageSource") {
                    cout << controller.addDCVoltageSource(name, node1, node2, value, circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                } else {
                    cout << controller.addDCCurrentSource(name, node1, node2, value, circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                }
            } else if (regex_match(input, match, addSINSource)) {
                if (match[1] != "V" && match[1] != "I") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number[3] = {match[5].str(), match[7].str(), match[9].str()};
                string unit[3] = {match[6].str(), match[8].str(), match[10].str()};
                double value[3];
                int check = 0;
                for (int i = 0; i < 3; i++) {
                    try {
                        value[i] = stod(number[i]);
                    } catch (const invalid_argument &e) {
                        cout << value[i] << endl;
                        cout << "Error: Invalid numeric value\n";
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                check = 0;
                for (int i = 0; i < 3; i++) {
                    string error = controller.handleError(name, circuit, value[i]);
                    if (error != "") {
                        cout << error;
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                for (int i = 0; i < 3; i++) {
                    if (!unit[i].empty()) {
                        switch (unit[i][0]) {
                            case 'G':
                                value[i] *= 1e9;
                                break;
                            case 'M':
                                value[i] *= 1e6;
                                break;
                            case 'k':
                            case 'K':
                                value[i] *= 1e3;
                                break;
                            case 'm':
                                value[i] *= 1e-3;
                                break;
                            case 'u':
                                value[i] *= 1e-6;
                                break;
                            case 'n':
                                value[i] *= 1e-9;
                                break;
                        }
                    }
                }
                if (value[2] <= 0) {
                    cout << "Frequency can not be Zero or Negative!\n";
                    continue;
                }
                if (match[1] == "V") {
                    cout << controller.addSinusoidalVoltageSource(name, node1, node2, value[1], value[2], value[0],
                                                                  circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                } else {
                    cout << controller.addSinusoidalCurrentSource(name, node1, node2, value[1], value[2], value[0],
                                                                  circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                }
            } else if (regex_match(input, match, addPULSESource)) {
                if (match[1] != "V" && match[1] != "I") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number[7] = {match[5].str(), match[7].str(), match[9].str(), match[11].str(), match[13].str(),
                                    match[15].str(), match[17].str()};
                string unit[7] = {match[6].str(), match[8].str(), match[10].str(),
                                  match[12].str(), match[14].str(), match[16].str(), match[18].str()};
                double value[7];
                int check = 0;
                for (int i = 0; i < 7; i++) {
                    try {
                        value[i] = stod(number[i]);
                    } catch (const invalid_argument &e) {
                        cout << value[i] << endl;
                        cout << "Error: Invalid numeric value\n";
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                check = 0;
                for (int i = 0; i < 7; i++) {
                    string error = controller.handleError(name, circuit, value[i]);
                    if (error != "") {
                        cout << error;
                        check = 1;
                        break;
                    }
                }
                if (check)
                    continue;
                for (int i = 0; i < 7; i++) {
                    if (!unit[i].empty()) {
                        switch (unit[i][0]) {
                            case 'G':
                                value[i] *= 1e9;
                                break;
                            case 'M':
                                value[i] *= 1e6;
                                break;
                            case 'k':
                            case 'K':
                                value[i] *= 1e3;
                                break;
                            case 'm':
                                value[i] *= 1e-3;
                                break;
                            case 'u':
                                value[i] *= 1e-6;
                                break;
                            case 'n':
                                value[i] *= 1e-9;
                                break;
                        }
                    }
                }
                if (value[2] < 0 || value[3] < 0 || value[4] < 0 || value[5] < 0 || value[6] < 0) {
                    cout << "Time parameters can not be Negative!\n";
                    continue;
                }
                if (value[6] < value[5] + value[4] + value[3] + value[2]) {
                    cout << "Period shorter than enough!\n";
                    continue;
                }
                if (match[1] == "V") {
                    cout << controller.addPulseVoltageSource(name, node1, node2, value[0], value[1], value[2], value[3],
                                                             value[4], value[5], value[6], circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                } else {
                    cout << controller.addPulseCurrentSource(name, node1, node2, value[0], value[1], value[2], value[3],
                                                             value[4], value[5], value[6], circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                }
            } else if (regex_match(input, match, add_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number = match[5].str();
                string unit = match[6].str();
                double value;
                try {
                    value = stod(number);
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                string error = controller.handleError(name, circuit, value);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (!unit.empty()) {
                    switch (unit[0]) {
                        case 'G':
                            value *= 1e9;
                            break;
                        case 'M':
                            value *= 1e6;
                            break;
                        case 'k':
                        case 'K':
                            value *= 1e3;
                            break;
                        case 'm':
                            value *= 1e-3;
                            break;
                        case 'u':
                            value *= 1e-6;
                            break;
                        case 'n':
                            value *= 1e-9;
                            break;
                    }
                }
                cout << controller.addNewElement(node1, node2, name, value, circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, remove_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C" && match[1] != "D" && match[1] != "V" && match[1] != "I" ) {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.removeElement(match[1].str() + match[2].str(), circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, add_diode)) {
                string name = match[1].str() + match[2].str();
                string node1 = match[3];
                string node2 = match[4];
                string model = match[5];
                string error = controller.handleError(name, circuit, 0);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (match[1] != "D") {
                    cout << "Error: Element " << name << " not found in library\n";
                    continue;
                }
                if (model == "D" || model == "Z") {
                    cout << controller.addDiode(node1, node2, name, model, circuit);
                    if (!loadedFilename.empty()) {
                        controller.saveCircuitToFile(circuit, loadedFilename);
                    }
                }
            } else if (regex_match(input, match, add_ground)) {
                if (match[1] != "GND") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.addGround(match[2], circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, delete_ground)) {
                if (match[1] != "GND") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.deleteGround(match[2], circuit);
            } else if (regex_match(input, match, show_nodes)) {
                controller.showNodes(circuit);
            } else if (regex_match(input, match, list)) {
                controller.list(circuit);
            } else if (regex_match(input, match, list_element)) {
                if (match[1] != "R" && match[1] != "C" && match[1] != "L" && match[1] != "D" && match[1] != "Z") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }
                controller.listElement(match[1], circuit);
            } else if (regex_match(input, match, rename_node)) {
                string old_name = match[1];
                string new_name = match[2];
                string err;
                circuit->renameNode(old_name, new_name, err);
                cout << err;
            } else if (regex_match(input, match, show_details)) {
                controller.showCircuitDetails(circuit);
            } else if (regex_match(input, match, new_file)) {
                string address = match[1];
                ifstream fin(address, ios::in);

                if (!fin) {
                    cout << "Error opening file!" << endl;
                    continue;
                }

                cout << "lets start" << endl;

                string line;
                regex word_regex("\\S+");
                while (getline(fin, line)) {
                    auto words_begin = sregex_iterator(line.begin(), line.end(), word_regex);
                    auto words_end = sregex_iterator();

                    vector<string> words;

                    for (auto it = words_begin; it != words_end; ++it) {
                        words.push_back(it->str());
                    }

                    string name;
                    string node1;
                    string node2;
                    string tmpValue;
                    string tmpAmplitude;
                    string tmpFrequency;
                    string tmpOffset;
                    string tmpVi;
                    string tmpVf;
                    string tmpTD;
                    string tmpTR;
                    string tmpTF;
                    string tmpTOn;
                    string tmpPeriod;


                    string type = words[0];
                    if (type == "R" || type == "L" || type == "C" || type == "I" || type == "V") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        tmpValue = words[4];
                        if (node1 == "GND") {
                            controller.addGround(node1, circuit);
                        }
                        if (node2 == "GND") {
                            controller.addGround(node2, circuit);
                        }

                        double value = 0;

                        if (tmpValue[0] == '-') {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpValue.back() == 'G') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                            } else if (tmpValue.back() == 'M') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                            } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                            } else if (tmpValue.back() == 'u') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                            } else if (tmpValue.back() == 'n') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                            } else if (tmpValue.back() == 'm') {
                                value = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                            } else {
                                value = stod(tmpValue);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (type == "R" || type == "C" || type == "L") {
                            cout << controller.addNewElement(node1, node2, name, value, circuit);

                        } else if (type == "V") {
                            cout << controller.addDCVoltageSource(name, node1, node2, value, circuit);
                        } else if (type == "I") {
                            cout << controller.addDCCurrentSource(name, node1, node2, value, circuit);
                        }
                    }
                    if (type == "D" || type == "Z") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        if (node1 == "GND") {
                            controller.addGround(node1, circuit);
                        }
                        if (node2 == "GND") {
                            controller.addGround(node2, circuit);
                        }
                        string model;
                        if (type == "D") {
                            model = "D";
                        }
                        if (type == "Z") {
                            model = "Z";
                        }
                        cout << controller.addDiode(node1, node2, name, model, circuit);
                    }
                    if (type == "VSIN" || type == "ISIN") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        tmpAmplitude = words[4];
                        tmpFrequency = words[5];
                        tmpOffset = words[6];
                        double amplitude;
                        double frequency;
                        double offset;
                        try {
                            if (tmpAmplitude.back() == 'G') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e9;
                            } else if (tmpAmplitude.back() == 'M') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e6;
                            } else if (tmpAmplitude.back() == 'k' || tmpAmplitude.back() == 'K') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e3;
                            } else if (tmpAmplitude.back() == 'u') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-6;
                            } else if (tmpAmplitude.back() == 'n') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-9;
                            } else if (tmpAmplitude.back() == 'm') {
                                amplitude = stod(tmpAmplitude.substr(0, tmpAmplitude.size() - 1)) * 1e-3;
                            } else {
                                amplitude = stod(tmpAmplitude);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpFrequency.back() == 'G') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e9;
                            } else if (tmpFrequency.back() == 'M') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e6;
                            } else if (tmpFrequency.back() == 'k' || tmpFrequency.back() == 'K') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e3;
                            } else if (tmpFrequency.back() == 'u') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-6;
                            } else if (tmpFrequency.back() == 'n') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-9;
                            } else if (tmpFrequency.back() == 'm') {
                                frequency = stod(tmpFrequency.substr(0, tmpFrequency.size() - 1)) * 1e-3;
                            } else {
                                frequency = stod(tmpFrequency);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        try {
                            if (tmpOffset.back() == 'G') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e9;
                            } else if (tmpOffset.back() == 'M') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e6;
                            } else if (tmpOffset.back() == 'k' || tmpOffset.back() == 'K') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e3;
                            } else if (tmpOffset.back() == 'u') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-6;
                            } else if (tmpOffset.back() == 'n') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-9;
                            } else if (tmpOffset.back() == 'm') {
                                offset = stod(tmpOffset.substr(0, tmpOffset.size() - 1)) * 1e-3;
                            } else {
                                offset = stod(tmpOffset);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        if (type == "VSIN") {
                            cout << controller.addSinusoidalVoltageSource(name, node1, node2, amplitude, frequency,
                                                                          offset, circuit);
                        }
                        if (type == "ISIN") {
                            cout << controller.addSinusoidalCurrentSource(name, node1, node2, amplitude, frequency,
                                                                          offset, circuit);
                        }
                    }

                    if (type == "VPULSE" || type == "IPULSE") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        tmpVi = words[4];
                        tmpVf = words[5];
                        tmpTD = words[6];
                        tmpTR = words[7];
                        tmpTF = words[8];
                        tmpTOn = words[9];
                        tmpPeriod = words[10];
                        double Vi, Vf, TD, TR, TF, TOn, period;
                        try {
                            if (tmpVi.back() == 'G') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e9;
                            } else if (tmpVi.back() == 'M') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e6;
                            } else if (tmpVi.back() == 'k' || tmpVi.back() == 'K') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e3;
                            } else if (tmpVi.back() == 'u') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-6;
                            } else if (tmpVi.back() == 'n') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-9;
                            } else if (tmpVi.back() == 'm') {
                                Vi = stod(tmpVi.substr(0, tmpVi.size() - 1)) * 1e-3;
                            } else {
                                Vi = stod(tmpVi);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpVf.back() == 'G') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e9;
                            } else if (tmpVf.back() == 'M') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e6;
                            } else if (tmpVf.back() == 'k' || tmpVf.back() == 'K') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e3;
                            } else if (tmpVf.back() == 'u') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-6;
                            } else if (tmpVf.back() == 'n') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-9;
                            } else if (tmpVf.back() == 'm') {
                                Vf = stod(tmpVf.substr(0, tmpVf.size() - 1)) * 1e-3;
                            } else {
                                Vf = stod(tmpVf);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTD.back() == 'G') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e9;
                            } else if (tmpTD.back() == 'M') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e6;
                            } else if (tmpTD.back() == 'k' || tmpTD.back() == 'K') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e3;
                            } else if (tmpTD.back() == 'u') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-6;
                            } else if (tmpTD.back() == 'n') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-9;
                            } else if (tmpTD.back() == 'm') {
                                TD = stod(tmpTD.substr(0, tmpTD.size() - 1)) * 1e-3;
                            } else {
                                TD = stod(tmpTD);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTR.back() == 'G') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e9;
                            } else if (tmpTR.back() == 'M') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e6;
                            } else if (tmpTR.back() == 'k' || tmpTR.back() == 'K') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e3;
                            } else if (tmpTR.back() == 'u') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-6;
                            } else if (tmpTR.back() == 'n') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-9;
                            } else if (tmpTR.back() == 'm') {
                                TR = stod(tmpTR.substr(0, tmpTR.size() - 1)) * 1e-3;
                            } else {
                                TR = stod(tmpTR);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTF.back() == 'G') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e9;
                            } else if (tmpTF.back() == 'M') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e6;
                            } else if (tmpTF.back() == 'k' || tmpTF.back() == 'K') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e3;
                            } else if (tmpTF.back() == 'u') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-6;
                            } else if (tmpTF.back() == 'n') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-9;
                            } else if (tmpTF.back() == 'm') {
                                TF = stod(tmpTF.substr(0, tmpTF.size() - 1)) * 1e-3;
                            } else {
                                TF = stod(tmpTF);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpTOn.back() == 'G') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e9;
                            } else if (tmpTOn.back() == 'M') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e6;
                            } else if (tmpTOn.back() == 'k' || tmpTOn.back() == 'K') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e3;
                            } else if (tmpTOn.back() == 'u') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-6;
                            } else if (tmpTOn.back() == 'n') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-9;
                            } else if (tmpTOn.back() == 'm') {
                                TOn = stod(tmpTOn.substr(0, tmpTOn.size() - 1)) * 1e-3;
                            } else {
                                TOn = stod(tmpTOn);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }

                        try {
                            if (tmpPeriod.back() == 'G') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e9;
                            } else if (tmpPeriod.back() == 'M') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e6;
                            } else if (tmpPeriod.back() == 'k' || tmpPeriod.back() == 'K') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e3;
                            } else if (tmpPeriod.back() == 'u') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-6;
                            } else if (tmpPeriod.back() == 'n') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-9;
                            } else if (tmpPeriod.back() == 'm') {
                                period = stod(tmpPeriod.substr(0, tmpPeriod.size() - 1)) * 1e-3;
                            } else {
                                period = stod(tmpPeriod);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (type == "VPULSE") {
                            cout << controller.addPulseVoltageSource(name, node1, node2, Vi, Vf, TD, TR, TF, TOn,
                                                                     period, circuit);
                        }
                        if (type == "IPULSE") {
                            cout << controller.addPulseCurrentSource(name, node1, node2, Vi, Vf, TD, TR, TF, TOn,
                                                                     period, circuit);
                        }
                    }
                    if (type == "E" || type == "G") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        string controlNode1 = words[4];
                        string controlNode2 = words[5];
                        string tmpGain = words[5];
                        double gain;
                        Node *cn1 = circuit->getCreateNode(controlNode1);
                        Node *cn2 = circuit->getCreateNode(controlNode2);
                        bool cn1Connected = false;
                        bool cn2Connected = false;
                        for (auto element: circuit->getElements()) {
                            if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                                cn1Connected = true;
                            }
                            if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                                cn2Connected = true;
                            }
                        }

                        if (!cn1Connected || !cn2Connected) {
                            cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected"
                                 << endl;
                            continue;
                        }
                        try {
                            if (tmpGain.back() == 'G') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e9;
                            } else if (tmpGain.back() == 'M') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e6;
                            } else if (tmpGain.back() == 'k' || tmpGain.back() == 'K') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e3;
                            } else if (tmpGain.back() == 'u') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-6;
                            } else if (tmpGain.back() == 'n') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-9;
                            } else if (tmpGain.back() == 'm') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-3;
                            } else {
                                gain = stod(tmpGain);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (gain == 0) {
                            cout << "Error: Gain cannot be zero" << endl;
                            continue;
                        }
                        if (type == "E") {
                            cout << controller.addVCVS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                        }
                        if (type == "G") {
                            cout << controller.addVCCS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                        }
                    }

                    if (type == "H" || type == "F") {
                        name = words[1];
                        node1 = words[2];
                        node2 = words[3];
                        string cElement = words[4];
                        string tmpGain = words[5];
                        double gain;

                        Element *controlElement = nullptr;
                        for (auto &element: circuit->getElements()) {
                            if (element->getName() == cElement) {
                                controlElement = element;
                                break;
                            }
                        }
                        if (controlElement == nullptr) {
                            cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                            continue;
                        }
                        try {
                            if (tmpGain.back() == 'G') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e9;
                            } else if (tmpGain.back() == 'M') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e6;
                            } else if (tmpGain.back() == 'k' || tmpGain.back() == 'K') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e3;
                            } else if (tmpGain.back() == 'u') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-6;
                            } else if (tmpGain.back() == 'n') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-9;
                            } else if (tmpGain.back() == 'm') {
                                gain = stod(tmpGain.substr(0, tmpGain.size() - 1)) * 1e-3;
                            } else {
                                gain = stod(tmpGain);
                            }
                        } catch (const invalid_argument &e) {
                            cout << "Error: Invalid numeric value for " << name << endl;
                            continue;
                        }
                        if (gain == 0) {
                            cout << "Error: Gain cannot be zero" << endl;
                            continue;
                        }
                        if (type == "H") {
                            cout << controller.addCCVS(name, node1, node2, cElement, gain, circuit);
                        }
                        if (type == "F") {
                            cout << controller.addCCCS(name, node1, node2, cElement, gain, circuit);
                        }

                    }

                }

                cout << "reading file ended :)" << endl;
                fin.close();
            } else if (regex_match(input, match, addVCVS)) {
                string name = 'E' + match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string controlNode1 = match[4].str();
                string controlNode2 = match[5].str();
                string tmpValue = match[6].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }
                Node *cn1 = circuit->getCreateNode(controlNode1);
                Node *cn2 = circuit->getCreateNode(controlNode2);
                bool cn1Connected = false;
                bool cn2Connected = false;
                for (auto element: circuit->getElements()) {
                    if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                        cn1Connected = true;
                    }
                    if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                        cn2Connected = true;
                    }
                }

                if (!cn1Connected || !cn2Connected) {
                    cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected";
                    continue;
                }
                cout << controller.addVCVS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, addVCCS)) {
                string name = 'G' + match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string controlNode1 = match[4].str();
                string controlNode2 = match[5].str();
                string tmpValue = match[6].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }
                Node *cn1 = circuit->getCreateNode(controlNode1);
                Node *cn2 = circuit->getCreateNode(controlNode2);
                bool cn1Connected = false;
                bool cn2Connected = false;
                for (auto element: circuit->getElements()) {
                    if (element->getFirstNode() == cn1 || element->getSecondNode() == cn1) {
                        cn1Connected = true;
                    }
                    if (element->getFirstNode() == cn2 || element->getSecondNode() == cn2) {
                        cn2Connected = true;
                    }
                }

                if (!cn1Connected || !cn2Connected) {
                    cout << "Error: " << controlNode1 << " and " << controlNode2 << " are not connected";
                    continue;
                }
                cout << controller.addVCCS(name, node1, node2, controlNode1, controlNode2, gain, circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, addCCVS)) {
                string name = 'H' + match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string cElement = match[4].str();
                string tmpValue = match[5].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }

                Element *controlElement = nullptr;
                for (auto &element: circuit->getElements()) {
                    if (element->getName() == cElement) {
                        controlElement = element;
                        break;
                    }
                }
                if (controlElement == nullptr) {
                    cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                    continue;
                }
                cout << controller.addCCVS(name, node1, node2, cElement, gain, circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, addCCCS)) {
                string name = 'F' + match[1].str();
                string node1 = match[2].str();
                string node2 = match[3].str();
                string cElement = match[4].str();
                string tmpValue = match[5].str();
                double gain;
                try {
                    if (tmpValue.back() == 'G') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e9;
                    } else if (tmpValue.back() == 'M') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e6;
                    } else if (tmpValue.back() == 'k' || tmpValue.back() == 'K') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e3;
                    } else if (tmpValue.back() == 'u') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-6;
                    } else if (tmpValue.back() == 'n') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-9;
                    } else if (tmpValue.back() == 'm') {
                        gain = stod(tmpValue.substr(0, tmpValue.size() - 1)) * 1e-3;
                    } else {
                        gain = stod(tmpValue);
                    }
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value for " << name << endl;
                    continue;
                }
                if (gain == 0) {
                    cout << "Error: Gain cannot be zero" << endl;
                    continue;
                }

                Element *controlElement = nullptr;
                for (auto &element: circuit->getElements()) {
                    if (element->getName() == cElement) {
                        controlElement = element;
                        break;
                    }
                }
                if (controlElement == nullptr) {
                    cout << "Error: " << cElement << " does not exist in the circuit" << endl;
                    continue;
                }
                cout << controller.addCCCS(name, node1, node2, cElement, gain, circuit);
                if (!loadedFilename.empty()) {
                    controller.saveCircuitToFile(circuit, loadedFilename);
                }
            } else if (regex_match(input, match, save_file)) {
                string name = match[1].str();
                name += ".txt";
                controller.saveCircuitToFile(circuit, name);
                cout << "Circuit saved to " << name << endl;
            } else if (regex_match(input, match, show_existing_schematics)) {
                controller.showAndLoadSchematic(circuit, loadedFilename);
            } else if (regex_match(input, match, preanalysischeck)) {
                controller.preAnalysisErrs(circuit);
            } else if (regex_match(input, match, check_matrix)) {
                controller.checkMatrix(circuit);
            } else if (regex_match(input, match, DC_Analaysis)) {
                if (!controller.preAnalysisErrs(circuit))
                    continue;
                controller.DCAnalysis(circuit);
            } else if (regex_match(input, match, DCSweep)) {
                string source = match[1];

                double start;
                double stop;
                double increament;
                string type = match[5].str();
                string nolement = match[6].str();
                try {
                    start = stod(match[2].str());
                    stop = stod(match[3].str());
                    increament = stod(match[4].str());
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                if (!controller.preAnalysisErrs(circuit))
                    continue;
                controller.DCSweep(source, start, stop, increament, type, nolement, circuit);
            }

            else if (regex_match(input, match, multipleDCSweep)){
                regex probePattern(R"([IV]\(\S+\))", regex::icase);
                sregex_iterator it(input.begin(), input.end(), probePattern);
                sregex_iterator end;
                vector<string> probes;
                while (it != end) {
                    probes.push_back(it->str());
                    ++it;
                }
                controller.multipleDCSweep(probes, match[1], stod(match[2]), stod(match[3]), stod(match[4]), circuit);
            }

            else if (regex_match(input, match, multipleTRAN)){
                regex probePattern(R"([IV]\(\S+\))", regex::icase);
                sregex_iterator it(input.begin(), input.end(), probePattern);
                sregex_iterator end;
                vector<string> probes;
                while (it != end) {
                    probes.push_back(it->str());
                    ++it;
                }
                controller.multipleTransient(probes, stod(match[1]), stod(match[2]), stod(match[3]), circuit);
            }

            else if (regex_match(input, match, transient)) {
                double start;
                double stop;
                double increament;
                string type = match[4].str();
                string nolement = match[5].str();
                try {
                    start = stod(match[1].str());
                    stop = stod(match[2].str());
                    increament = stod(match[3].str());
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                if (!controller.preAnalysisErrs(circuit))
                    continue;
                controller.TransientAnalysis(start, stop, increament, type, nolement, circuit);

            } else if (regex_match(input, match, reset)) {
                circuit->reset();
                cout << "A new schematic is created!\nStart again!\n";
            } else if (regex_match(input, match, exit)) {
                cout << "Bye Bye!\n";
                return;
            } else
                cout << "Syntax error\n";
        }
    }
};

int main() {
    View view;
    view.run();

    return 0;
}