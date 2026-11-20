#include "Glial.h"
#include "constants.h"
#include "voxel.h"

#include <cassert>
#include <cmath>
#include <vector>

namespace {

const double micrometers_to_internal = 1.0 / 1000.0;

Sphere makeSphere(int id, double x_um, double y_um, double z_um, double radius_um)
{
    return Sphere(id, 0,
                  Eigen::Vector3d(x_um, y_um, z_um) * micrometers_to_internal,
                  radius_um * micrometers_to_internal, 1);
}

void testPeriodicFaceEdgeAndCornerCrossings()
{
    const Eigen::Vector3d min_limits = Eigen::Vector3d::Zero();
    const Eigen::Vector3d max_limits = Eigen::Vector3d::Constant(0.1);

    Glial cell(0);
    std::vector<Sphere> geometry;
    geometry.push_back(makeSphere(0, 99.6009, 86.347, 34.1795, 0.90655));
    geometry.push_back(makeSphere(1, 99.6009, 99.6009, 34.1795, 0.90655));
    geometry.push_back(makeSphere(2, 99.6009, 99.6009, 99.6009, 0.90655));
    cell.set_spheres(geometry);
    cell.addPeriodicImages(min_limits, max_limits);

    Eigen::Vector3d face_point(0.001, 86.347, 34.1795);
    Eigen::Vector3d edge_point(0.001, 0.001, 34.1795);
    Eigen::Vector3d corner_point(0.001, 0.001, 0.001);
    face_point *= micrometers_to_internal;
    edge_point *= micrometers_to_internal;
    corner_point *= micrometers_to_internal;

    assert(cell.isPosInsideObstacle(face_point, 0.0));
    assert(cell.isPosInsideObstacle(edge_point, 0.0));
    assert(cell.isPosInsideObstacle(corner_point, 0.0));
    assert(cell.spheres.size() > geometry.size());
    assert(cell.spheres.back().object_id == 0);

    const Eigen::Vector3d start(99.999 * micrometers_to_internal,
                                86.347 * micrometers_to_internal,
                                34.1795 * micrometers_to_internal);
    Eigen::Vector3d step(1.0, 0.0, 0.0);
    const double full_step = 0.002 * micrometers_to_internal;
    Voxel voxel(min_limits, max_limits);
    Walker walker;
    walker.setVoxelPosition(start);
    walker.setRealPosition(start);
    walker.location = Walker::intra;
    walker.in_obj_type = glial_obstacle_type;
    walker.in_obj_index = 0;

    Collision membrane_before_wrap;
    assert(!cell.checkCollision(walker, step, full_step, membrane_before_wrap));
    assert(membrane_before_wrap.type == Collision::null);

    double face_step = full_step;
    Collision face_collision;
    assert(voxel.CheckCollision(walker, step, face_step, face_collision));
    assert(face_collision.type == Collision::hit);
    assert(face_collision.col_location == Collision::voxel);

    Eigen::Vector3d real_position = walker.pos_r;
    Eigen::Vector3d voxel_position = walker.pos_v;
    assert(mapPeriodicVoxelCrossing(real_position, voxel_position, step, face_collision.t,
                                    min_limits, max_limits));
    walker.setRealPosition(real_position);
    walker.setVoxelPosition(voxel_position);
    walker.status = Walker::on_voxel;

    double remaining_step = full_step - face_collision.t;
    Collision face_after_wrap;
    assert(!voxel.CheckCollision(walker, step, remaining_step, face_after_wrap));
    assert(face_after_wrap.type == Collision::null);

    Collision membrane_after_wrap;
    assert(!cell.checkCollision(walker, step, remaining_step, membrane_after_wrap));
    assert(membrane_after_wrap.type == Collision::null);
    assert(walker.location == Walker::intra && walker.in_obj_index == 0);
    assert(cell.isPosInsideObstacle(walker.pos_v, 0.0));

    walker.setRealPosition(walker.pos_r + remaining_step * step);
    walker.setVoxelPosition(walker.pos_v + remaining_step * step);
    const double unwrapped_displacement_um =
        (walker.pos_r.x() - start.x()) / micrometers_to_internal;
    assert(std::abs(unwrapped_displacement_um - 0.002) < 1e-9);
}

void testPeriodicImageKeepsMembranePermeability()
{
    const Eigen::Vector3d min_limits = Eigen::Vector3d::Zero();
    const Eigen::Vector3d max_limits = Eigen::Vector3d::Constant(0.1);
    Glial cell(0);
    std::vector<Sphere> geometry(1, makeSphere(0, 99.6009, 86.347, 34.1795, 0.90655));
    cell.set_spheres(geometry);
    cell.percolation = 1.0;
    cell.prob_cross_i_e = 1.0;
    cell.addPeriodicImages(min_limits, max_limits);

    Walker walker;
    walker.setVoxelPosition(Eigen::Vector3d(0.001, 86.347, 34.1795) * micrometers_to_internal);
    walker.location = Walker::intra;
    Eigen::Vector3d step(1.0, 0.0, 0.0);
    const double step_length = 0.7 * micrometers_to_internal;
    Collision collision;
    assert(cell.checkCollision(walker, step, step_length, collision));
    assert(collision.type == Collision::hit);
    assert(collision.perm_crossing > 0.0);
}

} // namespace

int main()
{
    testPeriodicFaceEdgeAndCornerCrossings();
    testPeriodicImageKeepsMembranePermeability();
}
