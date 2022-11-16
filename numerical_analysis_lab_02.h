#ifndef NUMERICAL_ANALYSIS_LAB_02_H
#define NUMERICAL_ANALYSIS_LAB_02_H

#include <QMainWindow>
#include "ui_numerical_analysis_lab_02.h"

QT_BEGIN_NAMESPACE
namespace Ui { class Numerical_analysis_lab_02; }
QT_END_NAMESPACE

class Numerical_analysis_lab_02 : public QMainWindow
{
    Q_OBJECT

public:
    Numerical_analysis_lab_02(QWidget *parent = Q_NULLPTR);
private slots:
    void plot_test_task();
    void plot_main_task();
private:
    Ui::Numerical_analysis_lab_02Class ui;
};

#endif // NUMERICAL_ANALYSIS_LAB_02_H
