
use crate::coord::datum;
use crate::coord::ellipsoid::Ellipsoid;
use crate::coord::error::ErrorCode;

pub struct WGS84 {
    ellipsoid: Ellipsoid,
}

impl datum::Datum for WGS84 {

    fn has_ellipsoid(&self) -> bool {
        true
    }

    fn ellipsoid(&self) -> Result<Ellipsoid, ErrorCode> {
        Ok( self.ellipsoid )
    }
}

impl WGS84 {

    pub fn as_epsg_4326() -> WGS84 {
        WGS84 {
            ellipsoid: Ellipsoid::wgs84(),
        }
    }
}