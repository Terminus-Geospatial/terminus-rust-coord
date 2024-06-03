
pub mod wgs84;

use crate::coord::ellipsoid::Ellipsoid;
use crate::coord::error::ErrorCode;


pub trait Datum {

    fn has_ellipsoid(&self) -> bool;

    fn ellipsoid(&self) -> Result<Ellipsoid,ErrorCode>;
}
