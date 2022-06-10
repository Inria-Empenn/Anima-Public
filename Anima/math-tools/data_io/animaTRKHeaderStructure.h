#pragma once

namespace anima
{

struct TRKHeaderStructure
{
    char id_string[6]; // ID string for track file. The first 5 characters must be "TRACK".
    short dim[3]; // Dimension of the image volume.
    float voxel_size[3]; // Voxel size of the image volume.
    float origin[3]; // Origin of the image volume. This field is not yet being used by TrackVis. That means the origin is always (0, 0, 0).
    short n_scalars; // Number of scalars saved at each track point (besides x, y and z coordinates).
    char scalar_name[10][20]; // Name of each scalar. Can not be longer than 20 characters each. Can only store up to 10 names.
    short n_properties; // Number of properties saved at each track.
    char property_name[10][20]; // Name of each property. Can not be longer than 20 characters each. Can only store up to 10 names.
    float vox_to_ras[4][4]; // matrix for voxel to RAS (crs to xyz) transformation. If vox_to_ras[3][3] is 0, it means the matrix is not recorded. This field is added from version 2.
    char reserved[444]; // Reserved space for future version.
    char voxel_order[4]; // Storing order of the original image data. Explained here.
    char pad2[4]; // Paddings.
    float image_orientation_patient[6]; // Image orientation of the original image. As defined in the DICOM header.
    char pad1[2]; // Paddings.
    unsigned char invert_x; // Inversion/rotation flags used to generate this track file. For internal use only.
    unsigned char invert_y; // As above.
    unsigned char invert_z; // As above.
    unsigned char swap_xy; // As above.
    unsigned char swap_yz; // As above.
    unsigned char swap_zx; // As above.
    int n_count; // Number of tracks stored in this track file. 0 means the number was NOT stored.
    int version; // Version number. Current version is 2.
    int hdr_size; // Size of the header. Used to determine byte swap. Should be 1000.
};

} // end namespace anima
