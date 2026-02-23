#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>
// точность
double Eps = 1e-8;

using namespace std;

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // бегущая сферическая волна
    double wave;
    // расстояние до источника
    double r;
    // амплитуды скоростей для лап
    double paw_vy;
    // Скорость
    double vx;
    double vy;
    double vz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), wave(0.0), r(0.0), paw_vy(0.0), vx(0.0), vy(0.0), vz(0.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double wave, double r, double paw_vy, double vx, double vy, double vz) 
            : x(x), y(y), z(z), wave(wave), r(r), paw_vy(paw_vy), vx(vx), vy(vy), vz(vz)
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};

// Класс элемента сетки
class Element
{
// Класс сетки будет friend-ом и элемента тоже
// (и вообще будет нагло считать его просто структурой)
friend class CalcMesh;

protected:
    // Индексы узлов, образующих этот элемент сетки
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> nodes;
    vector<Element> elements;

public:
    // Конструктор сетки из заданного stl-файла
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints, 
        const std::vector<double>& pawsAreas, double v0 = 0, 
        double k = 0, double omega = 0, double x0 = 0, double y0 = 0, double z0 = 0) {

        // Пройдём по узлам в модели gmsh
        nodes.resize(nodesCoords.size() / 3);
        for(unsigned int i = 0; i < nodesCoords.size() / 3; i++) {
            // Координаты заберём из gmsh
            double pointX = nodesCoords[i*3];
            double pointY = nodesCoords[i*3 + 1];
            double pointZ = nodesCoords[i*3 + 2];
            // положение плоской волны
            double r = pow(pow(pointX - x0, 2)+pow(pointY - y0, 2)+pow(pointZ - z0, 2), 0.5);
            double wave;
            if (r < Eps) { wave = 0; }
            else { wave = cos(k * r)/r; }

            // считаем амплитуды скоростей для лап
            double paw_vy = 0.0;
            for (unsigned int i = 0; i < pawsAreas.size() / 7; ++i) {
                if (pointX >= pawsAreas[7*i] && pointX <= pawsAreas[7*i+1] &&
                    pointY >= pawsAreas[7*i + 2] && pointY <= pawsAreas[7*i+3] &&
                    pointZ >= pawsAreas[7*i + 4] && pointZ <= pawsAreas[7*i+5]) {
                        if (i % 2) { paw_vy = 2e1*(pointZ - pawsAreas[7*i + 6]); } 
                        else {paw_vy = -2e1*(pointZ - pawsAreas[7*i + 6]); }
                    }
            }
            double vy = paw_vy + v0;
            nodes[i] = CalcNode(pointX, pointY, pointZ, wave, r, paw_vy, 0.0, vy, 0.0);
        }

        // Пройдём по элементам в модели gmsh
        elements.resize(tetrsPoints.size() / 4);
        for(unsigned int i = 0; i < tetrsPoints.size() / 4; i++) {
            elements[i].nodesIds[0] = tetrsPoints[i*4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i*4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i*4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i*4 + 3] - 1;
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, double step, double v0 = 0, double k = 0, double omega = 0) {
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // передвижение точек
            nodes[i].move(tau);
            // пробег волны
            if (nodes[i].r < Eps) { nodes[i].wave = 0; }
            else { nodes[i].wave = cos(k * nodes[i].r - omega * step * tau)/nodes[i].r; }
            // изменение скорости лапы
            nodes[i].vy = v0 + nodes[i].paw_vy * cos(omega * step * tau);
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto wave = vtkSmartPointer<vtkDoubleArray>::New();
        wave->SetName("wave");

        auto paw_velocity = vtkSmartPointer<vtkDoubleArray>::New();
        paw_velocity->SetName("paw_velocity");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        for(unsigned int i = 0; i < nodes.size(); i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(nodes[i].x, nodes[i].y, nodes[i].z);

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {nodes[i].vx, nodes[i].vy, nodes[i].vz};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            wave->InsertNextValue(nodes[i].wave);
            paw_velocity->InsertNextValue(nodes[i].paw_vy);
        }

        // Грузим точки в сетку
        unstructuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(wave);
        unstructuredGrid->GetPointData()->AddArray(paw_velocity);

        // А теперь пишем, как наши точки объединены в тетраэдры
        for(unsigned int i = 0; i < elements.size(); i++) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, elements[i].nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, elements[i].nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, elements[i].nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, elements[i].nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Создаём снапшот в файле с заданным именем
        string fileName = "../bug3d/bug3d_wave_move-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};

int main()
{
    // Шаг по времени
    double tau = 0.01;
    // количество итераций
    unsigned int stepnum = 100;
    // пространственная частота
    double k = 2 * M_PI / 25;
    // частота
    double omega = 2e1 / stepnum / tau;

    // расположение лап [x_min, x_max, y_min, y_max, z_min, z_max, z_attachment]
    std::vector<double> pawsAreas = {   75., 89., 138., 147., 1., 12., 7.,
                                        75., 87., 115., 124., 1., 10., 6.,
                                        76., 86., 94., 106., 3., 13., 6.,
                                        94., 108., 138., 147., 1., 12., 7.,
                                        98., 110., 115., 124., 1., 10., 6.,
                                        97., 107., 94., 106., 3., 13., 6.};

    // скорость жука
    double v0 = 1e2;

    const unsigned int GMSH_TETR_CODE = 4;

    // подгружаем готовую сетку из прошлой лабы
    gmsh::initialize();
    gmsh::open("../bug_low_mesh.msh");
    

    // Теперь извлечём из gmsh данные об узлах сетки
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);
    // ось жука
    double x0 = 0;
    for(unsigned int i = 0; i < nodesCoord.size() / 3; i++) {
        double x = nodesCoord[i*3];
        x0 = (x0 * i + x)/(i + 1);
    }
    double y0 = 170;
    double z0 = 35;
    
    // И данные об элементах сетки тоже извлечём, нам среди них нужны только тетраэдры, которыми залит объём
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // На всякий случай проверим, что номера узлов идут подряд и без пробелов
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Индексация в gmsh начинается с 1, а не с нуля. Ну штош, значит так.
        assert(i == nodeTags[i] - 1);
    }
    // И ещё проверим, что в тетраэдрах что-то похожее на правду лежит.
    assert(tetrsNodesTags->size() % 4 == 0);

    CalcMesh mesh(nodesCoord, *tetrsNodesTags, pawsAreas, v0, k, omega, x0, y0, z0);

    gmsh::finalize();

    mesh.snapshot(0);

    // Делаем шаги по времени, 
    // на каждом шаге считаем новое состояние и пишем его в VTK
    for(unsigned int step = 1; step < stepnum; step++) {
        mesh.doTimeStep(tau, step, v0, k, omega);
        mesh.snapshot(step);
    }
    return 0;
}
